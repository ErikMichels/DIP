cDIP <- function(new_data) {
  library(ggplot2)
  library(ggbeeswarm)
  library(reticulate)

  message("Please ensure TREM_1, IL_6, and Procalcitonin are in pg/ml, untransformed and unscaled.")
  message("This function uses Python. The first time use might take a few minutes. You might need to restart R afterwards.")

  # Step 1: Check for Python installation
  if (!py_available(initialize = FALSE)) {
    stop(
      "Python is not available. Please install Python from https://www.python.org/downloads/.\n",
      "During installation, ensure you check the option to 'Add Python to PATH'.\n",
      "After installation, restart R and try running this function again."
    )
  }

  message("Python is available. Proceeding with the function setup.")

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
    geom_beeswarm(cex = 3) +
    scale_color_gradient2(
      low = "#8BCCF1",
      mid = "#896DB0",
      high = "#006E78",
      midpoint = 0.5,
      breaks = c(0.0, 0.5, 1.0), # Include 0.0, 0.5, and 1.0 in the legend
      limits = c(0.0, 1.0),      # Ensure the scale includes these values
      name = "cDIP distribution in your cohort"              # Set legend title
    ) +
    theme_minimal() +
    coord_flip() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_blank(), # Remove y-axis title
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12), # Legend title style
      legend.text = element_text(size = 10), # Legend text style
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    )


  # Save results to the global environment
  assign("cDIP_results", results_df, envir = .GlobalEnv)

  return(results_df)
}
