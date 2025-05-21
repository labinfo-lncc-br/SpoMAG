#' Predict Sporulation Potential
#'
#' This function predicts the sporulation potential of MAGs using an ensemble learning model.
#' It uses probabilities from a Random Forest and SVM as inputs to a meta-model.
#'
#' @param binary_matrix A binary matrix (1/0) indicating gene presence/absence for each MAG. Must include a `genome_ID` column.
#'
#' @return A tibble with predicted class and probability of sporulation for each genome.
#' @import dplyr
#' @importFrom stats predict
#' @export
predict_sporulation <- function(binary_matrix) {
  # Locate model file from internal package data
  model_path <- system.file("extdata", "models_sporulation.RData", package = "SpoMAG")

  if (model_path == "" || !file.exists(model_path)) {
    stop("Model file not found. Please ensure 'models_sporulation.RData' is included in inst/extdata.")
  }

  # Load the trained models
  load(model_path)  # loads: rf_model, svm_model, meta_model

  required_packages <- c("caret", "kernlab", "randomForest")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required"))
    }
  }

  # Check if "genome_ID" exists
  if (!"genome_ID" %in% colnames(binary_matrix)) {
    stop("Input binary matrix must contain a 'genome_ID' column.")
  }

  # Remove identifiers
  features <- binary_matrix[, setdiff(names(binary_matrix), "genome_ID")]

  # Add missing predictors used in RF and SVM
  rf_vars <- setdiff(colnames(rf_model$trainingData), ".outcome")
  svm_vars <- setdiff(colnames(svm_model$trainingData), ".outcome")
  all_vars <- union(rf_vars, svm_vars)

  missing_vars <- setdiff(all_vars, colnames(features))
  for (var in missing_vars) {
    features[[var]] <- 0
  }

  # Reorder columns to match training
  features <- binary_matrix %>%
    dplyr::select(dplyr::all_of(union(rf_vars, svm_vars))) %>%
    tibble::as_tibble()

  # Generate predictions from base models
  prob_rf <- predict(rf_model, features, type = "prob")[, "Esporulante"]
  prob_svm <- predict(svm_model, features, type = "prob")[, "Esporulante"]

  # Create meta-model input
  meta_input <- data.frame(
    RF_Prob = prob_rf,
    SVM_Prob = prob_svm
  )

  # Final prediction using the meta-model
  prob_meta <- predict(meta_model, meta_input, type = "prob")[, "Esporulante"]
  class_meta <- ifelse(prob_meta > 0.5, "Sporulating", "Non_sporulating")

  # Combine results
  result <- dplyr::tibble(
    genome_ID = binary_matrix$genome_ID,
    RF_Prob = prob_rf,
    SVM_Prob = prob_svm,
    Meta_Prob_Sporulating = prob_meta,
    Meta_Prediction = class_meta
  )

  return(result)
}
