## ------------------------------------------------------------------------
library(PheCAP)

## ------------------------------------------------------------------------
dim(ehrdata)

## ------------------------------------------------------------------------
set.seed(1234)

## ------------------------------------------------------------------------
patient_id <- rownames(ehrdata)

## ------------------------------------------------------------------------
icd_feature <- data.frame(
  PatientID = patient_id,
  main_ICD = ehrdata[, "main_ICD"],
  ehrdata[, substr(colnames(ehrdata), 1, 3) == "COD"])

## ------------------------------------------------------------------------
nlp_feature <- data.frame(
  PatientID = patient_id,
  main_NLP = ehrdata[, "main_NLP"],
  ehrdata[, substr(colnames(ehrdata), 1, 3) == "NLP"])

## ------------------------------------------------------------------------
hu_feature <- data.frame(
  PatientID = patient_id,
  healthcare_utlization = ehrdata[, "healthcare_utlization"])

## ------------------------------------------------------------------------
id_label <- patient_id[!is.na(ehrdata$label)]
id_label_training <- sample(id_label, round(length(id_label)/2))
id_label_validation <- setdiff(id_label, id_label_training)

## ------------------------------------------------------------------------
training_label <- data.frame(
  PatientID = id_label_training, 
  label = ehrdata[patient_id %in% id_label_training, "label"])
table(training_label$label)

## ------------------------------------------------------------------------
validation_label <- data.frame(
  PatientID = id_label_validation, 
  label = ehrdata[patient_id %in% id_label_validation, "label"])
table(validation_label$label)

## ------------------------------------------------------------------------
data <- PhecapData(
  icd_feature = icd_feature,
  nlp_feature = nlp_feature,
  hu_feature = hu_feature,
  training_label = training_label,
  validation_label = validation_label,
  patient_index = "PatientID")

## ------------------------------------------------------------------------
surrogates <- list(
  PhecapSurrogate(variable_names = "main_ICD", 
                  lower_cutoff = 1, upper_cutoff = 10),
  PhecapSurrogate(variable_names = "main_NLP", 
                  lower_cutoff = 1, upper_cutoff = 10),
  PhecapSurrogate(variable_names = c("main_ICD", "main_NLP"), 
                  lower_cutoff = 1, upper_cutoff = 10))

## ---- cache=TRUE---------------------------------------------------------
feature_selected <- phecap_run_feature_extraction(data, surrogates)
print(feature_selected)

## ---- eval=FALSE---------------------------------------------------------
#  model <- phecap_train_phenotyping_model(data, surrogates, feature_selected)
#  # print(model)

## ---- eval=FALSE---------------------------------------------------------
#  validation <- phecap_validate_phenotyping_model(data, surrogates, model)
#  # print(validation)
#  phecap_plot_roc_curves(validation)

## ---- eval=FALSE---------------------------------------------------------
#  phenotype <- phecap_predict_phenotype(data, surrogates, model)

