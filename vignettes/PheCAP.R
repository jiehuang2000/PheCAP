## ------------------------------------------------------------------------
library(PheCAP)

## ---- eval=FALSE---------------------------------------------------------
#  data <- PhecapData(
#    icd_feature = "icd_feature.csv",
#    nlp_feature = "nlp_feature.csv",
#    hu_feature = "hu_feature.csv",
#    training_label = "training_label.csv",
#    validation_label = "validation_label.csv",
#    patient_index = "PatientID")

## ---- eval=FALSE---------------------------------------------------------
#  surrogates <- list(
#    PhecapSurrogate(
#      variable_names = "ICD1",
#      lower_cutoff = 1, upper_cutoff = 10),
#    PhecapSurrogate(
#      variable_names = "NLP1",
#      lower_cutoff = 1, upper_cutoff = 10),
#    PhecapSurrogate(
#      variable_names = c("ICD1", "NLP1"),
#      lower_cutoff = 1, upper_cutoff = 10))

## ---- eval=FALSE---------------------------------------------------------
#  feature_selected <- phecap_run_feature_extraction(
#    data, surrogates)
#  print(feature_selected)

## ---- eval=FALSE---------------------------------------------------------
#  model <- phecap_train_phenotyping_model(
#    data, surrogates, feature_selected)
#  print(model)

## ---- eval=FALSE---------------------------------------------------------
#  validation <- phecap_validate_phenotyping_model(
#    data, surrogates, model)
#  print(validation)
#  phecap_plot_roc_curves(validation)

## ---- eval=FALSE---------------------------------------------------------
#  phenotype <- phecap_predict_phenotype(
#    data, surrogates, model)

