## ------------------------------------------------------------------------
library(PheCAP)

## ------------------------------------------------------------------------
set.seed(123)
w <- rgamma(9000, 0.3)
icd_feature <- data.frame(
  PatientID = sample.int(10000, 9000),
  ICD1 = rpois(9000, 7 * (rgamma(9000, 0.2) + w) / 0.5),
  ICD2 = rpois(9000, 6 * (rgamma(9000, 0.8) + w) / 1.1),
  ICD3 = rpois(9000, 1 * rgamma(9000, 0.5) / 0.5),
  ICD4 = rpois(9000, 2 * rgamma(9000, 0.5) / 0.5))
w <- rgamma(9500, 0.4)
nlp_feature <- data.frame(
  PatientID = sample.int(10000, 9500),
  NLP1 = rpois(9500, 8 * (rgamma(9500, 0.2) + w) / 0.6),
  NLP2 = rpois(9500, 2 * (rgamma(9500, 1.1) + w) / 1.5),
  NLP3 = rpois(9500, 5 * (rgamma(9500, 0.1) + w) / 0.5),
  NLP4 = rpois(9500, 11 * rgamma(9500, 1.9) / 1.9),
  NLP5 = rpois(9500, 3 * rgamma(9500, 0.5) / 0.5),
  NLP6 = rpois(9500, 2 * rgamma(9500, 0.5) / 0.5),
  NLP7 = rpois(9500, 1 * rgamma(9500, 0.5) / 0.5))
hu_feature <- data.frame(
  PatientID = seq_len(10000),
  NoteCount = rpois(10000, 30 * rgamma(10000, 0.1) / 0.1))
df <- merge(
  merge(icd_feature, nlp_feature, by = "PatientID", all = TRUE),
  hu_feature, by = "PatientID", all = TRUE)
df[is.na(df)] <- 0
ii <- sample.int(5000, 400)
expr <- quote(plogis(
  -5 + 1.5 * log1p(ICD1) + log1p(NLP1) +
    0.8 * log1p(NLP4) - 0.5 * log1p(NoteCount)))
summary(with(df, eval(expr)))
training_label <- data.frame(
  PatientID = head(ii, 200),
  Label = rbinom(200, 1, with(df[head(ii, 200), ], eval(expr))))
validation_label <- data.frame(
  PatientID = tail(ii, 200),
  Label = rbinom(200, 1, with(df[tail(ii, 200), ], eval(expr))))

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
  PhecapSurrogate(variable_names = "ICD1", 
                  lower_cutoff = 1, upper_cutoff = 10),
  PhecapSurrogate(variable_names = "NLP1", 
                  lower_cutoff = 1, upper_cutoff = 10),
  PhecapSurrogate(variable_names = c("ICD1", "NLP1"), 
                  lower_cutoff = 1, upper_cutoff = 10))

## ------------------------------------------------------------------------
feature_selected <- phecap_run_feature_extraction(data, surrogates)
print(feature_selected)

## ------------------------------------------------------------------------
model <- phecap_train_phenotyping_model(data, surrogates, feature_selected)
# print(model)

## ------------------------------------------------------------------------
validation <- phecap_validate_phenotyping_model(data, surrogates, model)
# print(validation)
phecap_plot_roc_curves(validation)

## ------------------------------------------------------------------------
phenotype <- phecap_predict_phenotype(data, surrogates, model)

