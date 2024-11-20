cDIP <- function(new_data) {
  library(ggplot2)
  library(ggbeeswarm)
  library(reticulate)

  message("Please ensure TREM_1, IL_6, and Procalcitonin are in pg/ml, untransformed and unscaled.")
  message("This function uses Python. The first time use might take a few minutes. You might need to restart R afterwards.")

  # Define the base directory for the virtual environment
  package_base_dir <- system.file(package = "DIP")
  venv_dir <- file.path(package_base_dir, "r-reticulate-env")

  # Ensure the virtual environment exists
  if (!dir.exists(venv_dir)) {
    message("Creating a new virtual environment...")
    virtualenv_create(envname = venv_dir)
  }

  # Activate the virtual environment
  use_virtualenv(venv_dir, required = TRUE)

  # Define the required Python packages
  required_packages <- c("numpy", "pandas", "scikit-learn")

  # Install missing packages
  for (pkg in required_packages) {
    if (!py_module_available(pkg)) {
      message(sprintf("Installing Python package: %s", pkg))
      virtualenv_install(envname = venv_dir, packages = pkg)
    }
  }

  # Confirm dependencies are available
  if (!py_module_available("numpy") || !py_module_available("pandas") || !py_module_available("sklearn")) {
    stop("Required Python modules are missing even after installation. Check your Python setup.")
  }

  # Check if the model file exists
  model_path <- system.file("extdata/python", "model.pkl", package = "DIP")
  if (!file.exists(model_path)) {
    stop("Model file not found. Please check the package installation.")
  }

  # Load the model
  model <- py_load_object(model_path)

  # Validate input data
  if (!is.data.frame(new_data)) {
    stop("Error: new_data must be a data frame.")
  }
  if (!"ID" %in% names(new_data)) {
    stop("Error: Data must contain an 'ID' column.")
  }
  if (anyDuplicated(new_data$ID)) {
    stop("Error: IDs are not unique. Each ID must be unique.")
  }

  # Ensure required columns are present
  expected_vars <- c("Procalcitonin", "TREM_1", "IL_6")
  if (!all(expected_vars %in% names(new_data))) {
    stop("Error: Missing required columns: ", paste(setdiff(expected_vars, names(new_data)), collapse = ", "))
  }

  # Remove rows with missing classifier data
  rows_with_missing <- new_data$ID[rowSums(is.na(new_data[expected_vars])) > 0]
  if (length(rows_with_missing) > 0) {
    message("Patients with missing classifiers are omitted from the dataset. Affected patient IDs are:")
    message(paste(rows_with_missing, collapse = " "))
    new_data <- new_data[!new_data$ID %in% rows_with_missing, ]
  }

  # Prepare predictor data
  predictors <- new_data[, expected_vars]
  if (!all(sapply(predictors, is.numeric))) {
    stop("All predictor columns must be numeric.")
  }

  # Run the model prediction
  prediction <- model$model$predict(predictors)

  # Create results dataframe
  results_df <- data.frame(
    ID = new_data$ID,
    TREM_1 = new_data$TREM_1,
    IL_6 = new_data$IL_6,
    Procalcitonin = new_data$Procalcitonin,
    cDIP = prediction
  )

  # Plot the results
  p <- ggplot(results_df, aes(y = cDIP, x = factor(1), color = cDIP)) +
    geom_beeswarm(cex = 3, method = "hex", size = 1, dodge.width = 0.5) +
    scale_color_gradient2(low = "#8BCCF1", mid = "#896DB0", high = "#006E78", midpoint = 0.5) +
    theme_minimal() +
    scale_x_discrete(expand = expansion(add = 1.5)) +
    coord_flip() +
    ggtitle("cDIP Distribution in Your Cohort") +
    scale_color_gradient2(
      low = "#8BCCF1",
      mid = "#896DB0",
      high = "#006E78",
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    expand_limits(y = c(0, 1)) +
    theme(
      plot.title = element_text(size = 12, face = "italic"),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black")
    ) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"))
  suppressWarnings(print(p))

  # Save results to the global environment
  assign("cDIP_results", results_df, envir = .GlobalEnv)

  return(results_df)
}
