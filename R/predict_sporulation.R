#' Predict Sporulation Potential
#'
#' This function predicts the sporulation potential of MAGs using an ensemble learning model.
#' It uses probabilities from a Random Forest and SVM as inputs to a meta-model.
#'
#' @param binary_matrix A binary matrix (1/0) indicating gene presence/absence for each MAG. Must include a `codigo_genoma` column.
#' @param model_path Path to an .RData file containing trained objects: `rf_model`, `svm_model`, and `meta_model`.
#'
#' @return A tibble with RF and SVM probabilities, final prediction, and meta-model probability.
#' @export
#' @importFrom stats predict
#' @importFrom readr read_csv
utils::globalVariables(c(
  "%>%",
  "KEGG_ko", "Preferred_name", "codigo_genoma", "consensus_name_this_study",
  "presenca", "meta_model", "rf_model", "svm_model", "spoiiie_ko_eccca"
))

predict_sporulation <- function(binary_matrix, model_path) {
  # Load the trained models
  load(model_path)  # should load: rf_model, svm_model, meta_model

  required_packages <- c("caret", "kernlab", "randomForest")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required"))
    }
  }

  # Check if "codigo_genoma" exists
  if (!"codigo_genoma" %in% colnames(binary_matrix)) {
    stop("Input binary matrix must contain a 'codigo_genoma' column.")
  }

  # Remove identifiers
  features <- binary_matrix[, setdiff(names(binary_matrix), "codigo_genoma")]

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
  class_meta <- ifelse(prob_meta > 0.5, "Esporulante", "Nao_Esporulante")

  # Combine results
  result <- dplyr::tibble(
    codigo_genoma = binary_matrix$codigo_genoma,
    RF_Prob = prob_rf,
    SVM_Prob = prob_svm,
    Meta_Prediction = class_meta,
    Meta_Prob_Esporulante = prob_meta
  )

  return(result)
}
