# zzz.R

# Supress global variable notes in R CMD check
utils::globalVariables(c(
  "codigo_genoma", "consensus_name_this_study", "presenca",
  "rf_model", "svm_model", "meta_model",
  "Preferred_name", "KEGG_ko", "spoiiie_ko_eccca"
))
