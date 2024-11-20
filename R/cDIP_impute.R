#' Predicts the continuous degree of host response dysregulation: continuous Dysregulated Immune Profile. Imputes missing values.
#'
#' This function requires Python to be installed as it loads a complex trained Random Forest Regressor model.
#' The model predicts the degree of host response dysregulation on a continuous scale: the continuous Dysregulated Immune Profile or cDIP using plasma concentrations of sTREM-1, IL-6, and Procalcitonin.
#' The prediction model leverages absolute values to deliver precise, tailored outcomes, suitable for single-patient scenarios.
#' The function imputes missing values using 10 generated MICE datasets if necessary.
#' Be aware that this 3-biomarker based model with one imputed classifier has not been validated.
#' The continuous degree of host response dysregulation ranges from 0 (minor dysregulation) - 1 (major dysregulation).

#' @name cDIP_impute
#' @param new_data A data frame containing the following columns:
#' \itemize{
#'   \item \code{ID}: unique Identifier for the observations. Patients with multiple timepoints should have the timepoint included in their ID.
#'   \item \code{TREM_1}: sTREM-1 measured in picograms per milliliter (pg/ml). Data should be untransformed and unscaled. It is crucial that the
#'   column name is exactly \code{TREM_1} as the function's code explicitly references this name.
#'   \item \code{IL_6}: IL-6  measured in picograms per milliliter (pg/ml). Data should be untransformed and unscaled. It is crucial that the
#'   column name is exactly \code{IL_6} as the function's code explicitly references this name.
#'   \item \code{Procalcitonin}: Procalcitonin  measured in picograms per milliliter (pg/ml). Data should be untransformed and unscaled. It is crucial that the
#'   column name is exactly \code{Procalcitonin} as the function's code explicitly references this name.
#' }
#' @return A vector with predictions and a beeswarm plot
#' @import ggbeeswarm
#' @import ggplot2
#' @export
#' @examples
#' test_data <- data.frame(
#'   ID = 1:3,
#'   TREM_1 = c(182, 400, 1000),
#'   IL_6 = c(70, 5, 10000),
#'   Procalcitonin = c(877, 66, 20000)
#' )
#' result <- DIP_stage(test_data)
#' print(result)


cDIP_impute <- function(new_data) {
  library(ggplot2)
  library(ggbeeswarm)
  library(reticulate)

  message("Please ensure TREM_1, IL_6, and Procalcitonin are in pg/ml, untransformed and unscaled.")
  message("This function uses Python. The first time use might take a few minutes. You might need to restart R afterwards")

  # Define the base directory for the virtual environment, e.g., within the package installation directory
  package_base_dir <- system.file(package = "DIP")
  venv_dir <- file.path(package_base_dir, "r-reticulate-env")

  # Attempt to find a Python configuration automatically
  python_config <- py_discover_config()

  # If no Python installation is found, instruct the user to install Python
  if (is.null(python_config$python)) {
    stop("Error: No suitable Python installation found. This function requires Python. Please install Python at https://www.python.org/downloads/. During the installation make
         sure you tick the box of 'Add Python to PATH' prior to pressing 'Install now'. ")
  }

  # Create and configure a new virtual environment if it doesn't exist
  if (!dir.exists(venv_dir)) {
    message("Creating a new virtual environment for Python dependencies in", venv_dir, "...")
    virtualenv_create(envname = venv_dir, python = python_config$python)
  }

  # Activate the virtual environment
  use_virtualenv(venv_dir, required = TRUE)

  # Install required Python packages if they are missing
  required_packages <- c("numpy==1.24.4", "scikit-learn==1.4.2", "pandas==2.2.2")
  missing_packages <- required_packages[!sapply(required_packages, py_module_available)]
  if (length(missing_packages) > 0) {
    message("Installing required Python packages: ", paste(missing_packages, collapse=", "))
    virtualenv_install(envname = venv_dir, packages = missing_packages)
  }

  # Ensure that new_data is a data frame
  if (!is.data.frame(new_data)) {
    stop("Error: new_data must be a data frame.")
  }

  # Check for necessary data columns
  if (!"ID" %in% names(new_data)) {
    stop("Error: Data must contain an 'ID' column.")
  }
  if (anyDuplicated(new_data$ID)) {
    stop("Error: IDs are not unique. Each ID must be unique to ensure accurate data processing. Patients with multiple timepoints should have the timepoint included in their ID.")
  }

  # Locate and check the model file
  model_path <- system.file("extdata/python", "model.pkl", package = "DIP")
  if (!file.exists(model_path)) {
    stop("Model file not found. Please check the package installation.")
  }

  # Load the Python environment and model
  model <- py_load_object(model_path)

  # Ensure the predictor data is correct
  expected_vars <- c( "Procalcitonin", "TREM_1", "IL_6")
  if (!all(expected_vars %in% names(new_data))) {
    stop("Error: Not all required variables are present in the input data. Please use these exact columnnames for the respective biomarker
         'Procalcitonin', 'TREM_1', 'IL_6'")
  }

  # Run checks only if there are more than 10 patients
  if (nrow(new_data) > 10) {
    # Check for common value issues
    for (var in expected_vars) {
      # Calculate the frequency of each value
      freqs <- table(new_data[[var]])
      if (length(freqs) > 0) {
        # Find the maximum frequency
        max_freq <- max(freqs)
        percent_max_freq <- (max_freq / nrow(new_data)) * 100

        # Check if any value constitutes more than 10% of the data
        if (percent_max_freq > 10) {
          warning(sprintf("Warning: More than 10%% of %s are the exact same value (%.2f%%), suggesting a large number of samples are near the detection threshold.
          Numerous concentrations at the lower or upper limit can compromise the model’s performance.", var, percent_max_freq))
        }
      }
    }
  }

  # Check if there are more than 1 missing classifier
  rows_with_more_than_one_missing <- new_data$ID[rowSums(is.na(new_data[expected_vars])) > 1]
  if (length(rows_with_more_than_one_missing) > 1) {
    message("Patients missing more than 1 missing classifiers are omitted from the dataset. Affected patient IDs are:")
    message(paste(rows_with_more_than_one_missing, collapse = " "))
    new_data <- new_data[!new_data$ID %in% rows_with_more_than_one_missing, ]
  }

  # Calculate missing percentage per biomarker
  missing_percentages <- sapply(expected_vars, function(var) {
    sum(is.na(new_data[[var]])) / nrow(new_data) * 100
  })


  # Iterate over each biomarker and print missing percentage if applicable
  for (i in seq_along(expected_vars)) {
    var <- expected_vars[i]
    missing_percentage <- missing_percentages[i]

    if (missing_percentage > 0) {
      message(sprintf("Note: %.2f%% of %s values were missing and have been imputed using the median of 10 MICE imputed datasets", missing_percentage, var))
    }
  }


  imputed_count <- integer(nrow(new_data))  # To keep track of the number of imputations

  for (var in expected_vars) {
    missing_percentage <- sum(is.na(new_data[[var]])) / nrow(new_data) * 100
    if (missing_percentage > 0) {
      if (missing_percentage > 20) {
        warning(sprintf("Warning: Over 20%% of %s data is missing (%.2f%%) imputed values might not reflect the underlying distribution.", var, missing_percentage))
      }

      missing_indices <- which(is.na(new_data[[var]]))
      # Update the imputed count for rows where the variable was imputed
      imputed_count[missing_indices] <- imputed_count[missing_indices] + 1

      temp_data <- new_data[ , c("ID", "TREM_1", "IL_6", "Procalcitonin"), drop = FALSE]
      imp <- mice(temp_data[ , -1], pri=FALSE, maxit=0)
      predM <- imp$predictorMatrix
      imputed_sets <- mice(data=temp_data, m=10, maxit=10, meth="cart", pred=predM,  print = F)

      # Get completed data
      imp_long <- complete(imputed_sets, action = "long", include = FALSE)
      median_imputed <- aggregate(imp_long[, var] ~ ID, imp_long, median)

      # Iterate over each biomarker and calculate the median value for imputed biomarker across all imputations for each patient
      for (i in 1:length( expected_vars )) {
        # Calculate the median value for imputed biomarker across all imputations for each patient
        median_imputed <- aggregate(imp_long[,  expected_vars [i]] ~ ID, imp_long, median)
        # Replace missing values with median imputed values
        new_data <- new_data[order(match(new_data$ID, median_imputed$ID)), ]
        new_data[ expected_vars [i]] <- median_imputed[, 2]
      }
    }

    # Create labels for the "Warning" column based on the count of biomarkers imputed
    warning_labels <- rowSums(is.na(new_data[expected_vars]))
    warning_labels <- ifelse(imputed_count == 0, "NA",
                             ifelse(imputed_count == 1, "1 classifier imputed","Error"))

  }

  predictors <- new_data[, expected_vars, drop = FALSE]
  if (!all(sapply(predictors, is.numeric))) {
    stop("All predictor columns in new_data must be numeric.")
  }

  # Extract ID column and predictor columns, ensuring order
  ids <- new_data$ID
  predictors <- new_data[, expected_vars]

  # Convert factors to numeric codes if present
  predictors[] <- lapply(predictors, function(x) if(is.factor(x)) as.numeric(as.factor(x)) else x)

  # Ensure all predictor data is numeric
  if (!all(sapply(predictors, is.numeric))) {
    stop("All predictor columns in new_data must be numeric.")
  }

  # Prepare data for prediction
  predictors_matrix <- predictors
  colnames(predictors_matrix) <- expected_vars

  # Make and return the prediction
  prediction <- model$model$predict(predictors_matrix)

  # Create a dataframe for results
  results_df <- data.frame(
    ID = new_data$ID,
    TREM_1 = new_data$TREM_1,
    IL_6 = new_data$IL_6,
    Procalcitonin = new_data$Procalcitonin,
    cDIP = prediction,
    Warning = warning_labels)

  # Create beeswarm plot
  p <- ggplot(results_df, aes(y = cDIP, x = factor(1), color = cDIP)) +
    geom_beeswarm(cex = 3, method = "hex", groupOnX = FALSE, size = 1,   dodge.width = 0.5) + # Adjust method and groupOnX
    scale_color_gradient2(low = "#8BCCF1", mid = "#896DB0", high = "#006E78", midpoint = 0.5) +
    theme_minimal() +
    scale_x_discrete(expand = expansion(add = 1.5)) + # Add padding around factor levels
    coord_flip() +
    ggtitle("cDIP distrtibution in your Cohort") +
    scale_color_gradient2(
      low = "#8BCCF1",
      mid = "#896DB0",
      high = "#006E78",
      midpoint = 0.5,
      limits = c(0, 1)  # Force legend to range from 0 to 1

    ) +
    expand_limits(y = c(0, 1)) +  # Ensure all dots fit within 0 to 1
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
  print(p)


  # Save the dataframe to the environment
  assign("cDIP_impute_results", results_df, envir = .GlobalEnv)


  detach("package:reticulate", unload=TRUE)

  return(results_df)

}



