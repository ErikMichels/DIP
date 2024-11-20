#' Predicts the stage of host response dysregulation: Dysregulated Immune Profile 1-3. Imputes missing values.
#'
#' This function uses an extreme XGBoost model to predict the degree of host response dysregulation using using plasma concentrations of sTREM-1, IL-6, and Procalcitonin.
#' The prediction model leverages absolute values to deliver precise, tailored outcomes, suitable for single-patient scenarios.
#' The function imputes missing values using 10 generated MICE datasets if necessary.
#' Be aware that this 3-biomarker based model with one imputed classifier has not been validated.
#' Outputs include the original classifier data, the predicted DIP stage, and probability scores for each predicted stage.
#' DIP1 represent minor, DIP2 moderate and DIP3 major host response dysregulation.
#' Interactive visualizations are also provided when run in an appropriate environment.
#'
#'
#' @name DIP_stage_impute
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
#' @return A data frame with the original data alongside predictions and probability scores for each stage.
#' @import xgboost
#' @import plotly
#' @import ggplot2
#' @import mice
#' @importFrom plotly plot_ly
#' @importFrom plotly last_plot
#' @importFrom stats predict
#' @importFrom utils install.packages
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


DIP_stage_impute <- function(new_data) {
  ## Ensure required packages are loaded
  library(xgboost)
  library(plotly)
  library(ggplot2)
  library(mice)

  # Load the model from the package's internal directory
  model_path <- system.file("extdata", "xgb_model.json", package = "DIP")
  if (!file.exists(model_path)) {
    stop("Model file not found in package.")
  }
  model <- xgboost::xgb.load(model_path)

  message("Please ensure TREM_1, IL_6, and Procalcitonin are in pg/ml, untransformed and unscaled.")

  # Ensure that new_data is a data frame
  if (!is.data.frame(new_data)) {
    stop("Error: new_data must be a data frame.")
  }

  # Check if 'ID' column is present
  if (!"ID" %in% names(new_data)) {
    stop("Error: Data must contain an 'ID' column.")
  }

  # Check for unique IDs
  if (anyDuplicated(new_data$ID) > 0) {
    stop("Error: IDs are not unique. Each ID must be unique to ensure accurate data processing. Patients with multiple timepoints should have
         the timepoint included in their ID")
  }

  # Ensure the predictor data is correct
  expected_vars <- c( "Procalcitonin", "TREM_1", "IL_6")
  if (!all(expected_vars %in% names(new_data))) {
    stop("Error: Not all required variables are present in the input data. Please use these exact column names for the respective biomarker
         'Procalcitonin', 'TREM_1', 'IL_6'")
  }

  # Define the expected order of the predictor columns
  expected_vars <- c("TREM_1", "IL_6", "Procalcitonin")

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

  # Calculate missing percentage per biomarker
  missing_percentages <- sapply(expected_vars, function(var) {
    sum(is.na(new_data[[var]])) / nrow(new_data) * 100
  })

  # Check if there are more than 1 missing classifier
  rows_with_more_than_one_missing <- new_data$ID[rowSums(is.na(new_data[expected_vars])) > 1]
  if (length(rows_with_more_than_one_missing) > 1) {
    message("Patients missing more than 1 missing classifiers are omitted from the dataset. Affected patient IDs are:")
    message(paste(rows_with_more_than_one_missing, collapse = " "))
    new_data <- new_data[!new_data$ID %in% rows_with_more_than_one_missing, ]
  }


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

  # Extract ID column and predictor columns, ensuring order
  ids <- new_data$ID
  predictors <- new_data[, expected_vars]

  # Convert factors to numeric codes if present
  predictors[] <- lapply(predictors, function(x) if(is.factor(x)) as.numeric(as.factor(x)) else x)

  # Ensure all predictor data is numeric
  if (!all(sapply(predictors, is.numeric))) {
    stop("All predictor columns in new_data must be numeric.")
  }

  # Convert to matrix, as xgb.DMatrix requires a matrix or a numeric vector
  predictors_matrix <- as.matrix(predictors)

  # Create DMatrix object
  dmatrix <- xgboost::xgb.DMatrix(predictors_matrix)

  # Predict using the loaded model
  predictions <- predict(model, dmatrix, predcontrib = FALSE)

  # Convert probabilities to predicted classes
  num_classes <- 3  # This must match the number of classes your model predicts
  predicted_classes <- matrix(predictions, ncol = num_classes, byrow = TRUE)

  # Find the class index with the maximum probability
  max_indices <- apply(predicted_classes, 1, which.max)

  # Map numerical predictions to class labels
  labels <- c("DIP1", "DIP2", "DIP3")
  labeled_predictions <- labels[max_indices]

  # Create a dataframe for results
  results_df <- data.frame(
    ID = ids,
    TREM_1 = new_data$TREM_1,
    IL_6 = new_data$IL_6,
    Procalcitonin = new_data$Procalcitonin,
    DIP = labeled_predictions,
    DIP1_Prob = predicted_classes[,1],
    DIP2_Prob = predicted_classes[,2],
    DIP3_Prob = predicted_classes[,3],
    Warning = warning_labels  # Number of biomarkers imputed
  )

  # Create 3D scatter plot using plotly
  scatter_3d <- plot_ly(data = results_df, x = ~DIP1_Prob * 100, y = ~DIP2_Prob * 100, z = ~DIP3_Prob * 100, type = 'scatter3d', mode = 'markers',
                        color = ~DIP,
                        colors = c("#8BCCF1", "#896DB0", "#006E78"),
                        size = 2,
                        sizes = 100,
                        alpha = 1) %>%
    layout(scene = list(xaxis = list(title = "DIP1 Prob (%)"),
                        yaxis = list(title = "DIP2 Prob (%)"),
                        zaxis = list(title = "DIP3 Prob (%)")))
  print(scatter_3d)  # Display the 3D scatter plot


  # Wait for user to press Enter to continue
  cat("Press Enter to exchange the interactive 3D scatterplot for a piechart")
  readline()

  ## make piechart
  prediction_counts <- table(results_df$DIP)  # Count the frequency of each prediction

  # Convert to data frame and calculate percentages
  prediction_data <- as.data.frame(prediction_counts)
  prediction_data$percentage <- round(prediction_data$Freq / sum(prediction_data$Freq) * 100, 1)  # Calculate percentage
  prediction_data$fill_label <- paste(prediction_data$Var1, paste0("(", prediction_data$percentage, "%", ")"))  # Create label for plot and legend

  # Plotting using ggplot2
  colors <-  c("#8BCCF1", "#896DB0", "#006E78")
  # Plotting pie chart using ggplot2 and converting to plotly
  pie_chart <- ggplot2::ggplot(prediction_data, aes(x = "", y = Freq, fill = fill_label)) +
    ggplot2::geom_bar(width = 1, stat = "identity") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors, labels = prediction_data$fill_label) +
    theme_void() +
    labs(fill = "Prediction", title = "Distribution of DIP Predictions")
  print(pie_chart)  # Display the pie chart

  # Save the dataframe to the environment
  assign("DIP_stage_imputed", results_df, envir = .GlobalEnv)

  # Return results dataframe as the function's primary output
  return(results_df)
}
