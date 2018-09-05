phecap_read_or_set_frame <- function(source)
{
  data <- NULL
  if (is.character(source)) {
    if (endsWith(source, ".csv")) {
      data <- read.csv(
        source, header = TRUE, sep = ",",
        stringsAsFactors = FALSE)
    } else {
      data <- read.table(
        source, header = TRUE, sep = "\t",
        stringsAsFactors = FALSE)
    }
  } else {
    data <- as.data.frame(source)
  }
  if (nrow(data) == 0L || ncol(data) == 0L) {
    data <- NULL
  }
  
  return(data)
}


#' Define or Read Datasets for Phenotyping
#' 
#' Specify the data to be used for phenotyping.
#' 
#' Note that icd_feature, nlp_feature and hu_feature cannot be NULL, while training_label
#' can be NULL if only feature extraction is needed, and validation_label can be NULL 
#' if only feature extraction and model training is needed.
#' 
#' @param icd_feature a character scalar of the path of ICD feature, 
#' or a data.frame of the actual ICD feature.
#' @param nlp_feature a character scalar of the path of NLP feature, 
#' or a data.frame of the actual NLP feature.
#' @param hu_feature a character scalar of the path of HU feature, 
#' or a data.frame of the actual HU feature. One may also include other 
#' variables that must be included in the model.
#' @param training_label a character scalar of the path of training label, 
#' or a data.frame of the actual training label, or NULL.
#' @param validation_label a character scalar of the path of validation label, 
#' or a data.frame of the actual validation label, or NULL.
#' @param patient_index a character scalar or vector that specify a unique identifier 
#' for each patient in the data. Such variable(s) must appear in all the dataset above.
#' @param feature_transformation a function that will be applied to all the features. 
#' Since count data are typically right-skewed, by default \code{log1p} will be used.
#' 
#' @return An object of class \code{PhecapData}.
#' 
#' @export
PhecapData <- function(
  icd_feature, nlp_feature, hu_feature,
  training_label = NULL, validation_label = NULL,
  patient_index = "PatientID",
  feature_transformation = log1p)
{
  icd_feature <- phecap_read_or_set_frame(icd_feature)
  nlp_feature <- phecap_read_or_set_frame(nlp_feature)
  hu_feature <- phecap_read_or_set_frame(hu_feature)
  if (is.null(icd_feature)) {
    stop("ICD feature not found")
  }
  if (is.null(nlp_feature)) {
    stop("NLP feature not found")
  }
  if (is.null(hu_feature)) {
    stop("HU feature not found")
  }
  
  training_label <- phecap_read_or_set_frame(training_label)
  validation_label <- phecap_read_or_set_frame(validation_label)
  
  if (!(patient_index %in% names(icd_feature))) {
    stop(sprintf("ICD feature doesn't contain %s", patient_index))
  }
  if (!(patient_index %in% names(nlp_feature))) {
    stop(sprintf("NLP feature doesn't contain %s", patient_index))
  }
  if (!(patient_index %in% names(hu_feature))) {
    stop(sprintf("HU feature doesn't contain %s", patient_index))
  }

  if (!is.null(training_label) && 
      !(patient_index %in% names(training_label))) {
    stop(sprintf("Training label doesn't contain %s", patient_index))
  }
  if (!is.null(validation_label) && 
      !(patient_index %in% names(validation_label))) {
    stop(sprintf("Validation label doesn't contain %s", patient_index))
  }
  
  if (!is.function(feature_transformation)) {
    stop("feature_transformation should be a function")
  }
  
  data <- list(
    icd_feature = icd_feature,
    nlp_feature = nlp_feature,
    hu_feature = hu_feature,
    training_label = training_label,
    validation_label = validation_label,
    patient_index = patient_index,
    feature_transformation = feature_transformation)
  class(data) <- "PhecapData"
  
  return(data)
}

#' Define a Surrogate Variable used in Surrogate-Assisted Feature Extraction (SAFE)
#'
#' Define a surrogate varible from existing features, and specify associated lower and upper cutoffs.
#' 
#' This function only stores the definition. No calculation is done.
#' 
#' @usage PhecapSurrogate(variable_names, lower_cutoff = 1L, upper_cutoff = 10L)
#' 
#' @param variable_names a character scalar or vector consisting of variable names. If a vector is given, 
#' the value of the surrogate is defined as the sum of the values of each variable.
#' @param lower_cutoff a numeric scalar. If the surrogate value of a patient is less than or equal to 
#' this cutoff, then this patient is treated as a control in SAFE.
#' @param upper_cutoff a numeric scalar. If the surrogate value of a patient is greater than or equal to 
#' this cutoff, then this patient is treated as a case in SAFE.
#' 
#' @return An object of class \code{PhecapSurrogate}.
#' 
#' @export
PhecapSurrogate <- function(
  variable_names, 
  lower_cutoff = 1L, upper_cutoff = 10L)
{
  result <- list(
    variable_names = variable_names,
    lower_cutoff = lower_cutoff,
    upper_cutoff = upper_cutoff)
  class(result) <- "PhecapSurrogate"
  
  return(result)
}


phecap_check_surrogates <- function(
  surrogates, variable_list)
{
  for (surrogate in surrogates) {
    if (!all(surrogate$variable_names %in% variable_list)) {
      stop(sprintf("Variable(s) %s not found", paste0(
        setdiff(surrogate$variable_names, variable_list), 
        collapse = ", ")))
    }
  }
}


#' Run Surrogate-Assisted Feature Extraction (SAFE)
#' 
#' Run surrogate-assisted feature extraction (SAFE) using unlabeled data and subsampling.
#' 
#' In this unlabeled setting, the extremes of each surrogate are used to define cases and controls.
#' The variables selected are those selected in at least half of the subsamples.
#' 
#' @usage phecap_run_feature_extraction(data, surrogates,
#' subsample_size = 1000L, num_subsamples = 200L,
#' start_seed = 123L, verbose = 10L)
#' 
#' @param data An object of class PhecapData, obtained by calling PhecapData(...)
#' @param surrogates A list of objects of class PhecapSurrogate, obtained by something like
#' list(PhecapSurrogate(...), PhecapSurrogate(...))
#' @param subsample_size An integer scalar giving the size of each subsample
#' @param num_subsamples The number of subsamples drawn for each surrogate
#' @param start_seed in the i-th subsample, the seed is set to start_seed + i.
#' @param verbose print progress every \code{verbose} subsample if \code{verbose} is positive,
#' or remain quiet if \code{verbose} is zero
#' 
#' @return A character vector consisting of the names of the variables selected,
#' with attribute frequency, which consists of the proportion of being selected
#' among all subsamples for each variable.
#' 
#' @export
phecap_run_feature_extraction <- function(
  data, surrogates, 
  subsample_size = 1000L, num_subsamples = 200L, 
  start_seed = 123L, verbose = 10L)
{
  frame <- merge(
    data$icd_feature, data$nlp_feature, 
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$hu_feature,
    by = data$patient_index, all = TRUE)
  
  variable_list <- setdiff(
    names(frame), data$patient_index)
  variable_matrix <- as.matrix(frame[, variable_list])
  variable_matrix[is.na(variable_matrix)] <- 0.0
  phecap_check_surrogates(surrogates, variable_list)

  extremes <- lapply(
    surrogates, function(surrogate) {
      x <- rowSums(variable_matrix[, surrogate$variable_names, 
                                   drop = FALSE])
      cases <- which(x >= surrogate$upper_cutoff)
      controls  <- which(x <= surrogate$lower_cutoff)
      if (length(cases) <= subsample_size %/% 2L) {
        stop(sprintf("'%s' has too few cases; %s or %s",
             paste0(surrogate$variable_names, collapse = "&"),
             "decrease upper_cutoff", "decrease subsample_size"))
      }
      if (length(controls) <= subsample_size %/% 2L) {
        stop(sprintf("'%s' has too few controls; %s or %s", 
             paste0(surrogate$variable_names, collapse = "&"),
             "increase lower_cutoff", "decrease subsample_size"))
      }
      list(cases = cases, controls = controls)
    })
  if (verbose > 0L) {
    message <- sapply(
      extremes, function(z) 
        c(NumCases = length(z$cases),
          NumControls = length(z$controls)))
    colnames(message) <- sapply(surrogates, function(surrogate)
      paste0(surrogate$variable_names, collapse = "&"))
    print(message)
  }

  variable_matrix <- data$feature_transformation(variable_matrix)
  selection <- lapply(seq_along(surrogates), function(k) {
    surrogate <- surrogates[[k]]
    exclusion <- match(surrogate$variable_names, variable_list)
    cases <- extremes[[k]]$cases
    controls <- extremes[[k]]$controls
    half_size <- subsample_size %/% 2L
    if (verbose > 0L) {
      cat("Using surrogate", 
          paste0(surrogate$variable_names, collapse = "&"), "\n")
    }
    
    # subsampling
    nonzero <- t(sapply(seq_len(num_subsamples), function(ss) {
      if (verbose > 0L && (ss %% verbose == 1L || ss == num_subsamples)) {
        cat("Subsample", ss, "\n")
      }
      set.seed(start_seed + ss)
      ipos <- sort.int(sample(cases, subsample_size %/% 2L, FALSE))
      ineg <- sort.int(sample(controls, subsample_size %/% 2L, FALSE))
      
      y <- c(rep(1.0, length(ipos)), rep(0.0, length(ineg)))
      x <- variable_matrix[c(ipos, ineg), -exclusion, drop = FALSE]
      
      alpha <- fit_lasso_bic(x, y)
      alpha <- alpha[-1L]  # drop intercept
      as.integer(alpha != 0.0)
    }))
    
    nonzero_final <- matrix(1, nrow(nonzero), ncol(variable_matrix))
    nonzero_final[, -exclusion] <- nonzero
    nonzero_final
  })
  
  selection <- do.call("rbind", selection)
  frequency <- colMeans(selection)
  names(frequency) <- variable_list
  selected <- variable_list[frequency >= 0.5]

  result <- selected
  attr(result, "frequency") <- frequency
  class(result) <- "PhecapFeatureExtraction"
  result
}


phecap_generate_feature_matrix <- function(
  data, surrogates, frame, feature_selected)
{
  surrogate_matrix <- sapply(surrogates, function(surrogate) {
    rowSums(frame[, surrogate$variable_names, drop = FALSE])
  })
  surrogate_matrix <- data$feature_transformation(surrogate_matrix)
  colnames(surrogate_matrix) <- sapply(surrogates, function(surrogate) {
    paste0(surrogate$variable_names, collapse = "&")
  })
  surrogate_matrix[is.na(surrogate_matrix)] <- 0.0
  
  hu_matrix <- as.matrix(
    frame[, setdiff(
      names(data$hu_feature), data$patient_index),
      drop = FALSE])
  hu_matrix <- data$feature_transformation(hu_matrix)
  hu_matrix[is.na(hu_matrix)] <- 0.0
  
  other_matrix <- as.matrix(
    frame[, setdiff(feature_selected, c(
      colnames(surrogate_matrix), colnames(hu_matrix))),
      drop = FALSE])
  other_matrix <- data$feature_transformation(other_matrix)
  other_matrix[is.na(other_matrix)] <- 0.0
  other_matrix <- qr.resid(qr(cbind(
    1.0, surrogate_matrix, hu_matrix)), other_matrix)
  
  result <- cbind(
    surrogate_matrix,
    hu_matrix,
    other_matrix)
  attr(result, "free") <- ncol(surrogate_matrix) + ncol(hu_matrix)
  return(result)
}


#' Train Phenotyping Model using the Training Labels
#' 
#' Train the phenotyping model on the training dataset, and evaluate 
#' its performance via random splits of the training dataset.
#' 
#' @usage phecap_train_phenotyping_model(
#' data, surrogates, feature_selected,
#' method = "lasso_bic", train_percent = 0.7,
#' num_splits = 200L, start_seed = 12345L, verbose = 10L)
#' 
#' @param data an object of class \code{PhecapData}, obtained by 
#' calling \code{PhecapData(...)}.
#' 
#' @param surrogates a list of objects of class \code{PhecapSurrogate}, obtained 
#' by something like \code{list(PhecapSurrogate(...), PhecapSurrogate(...))}.
#' 
#' @param feature_selected a character vector of the features that should be included 
#' in the model, probably returned by \code{phecap_run_feature_extraction}.
#' 
#' @param method Either a character scalar or a list of two components.
#' If a character scalar is used, possible values are
#' \code{'plain'} (logistic regression without penalty),
#' \code{'lasso_cv'} (logistic regression with lasso penalty and CV tuning),
#' \code{'lasso_bic'} (logistic regression with lasso penalty and BIC tuning),
#' \code{'svm'} (support vector machine with CV tuning, package \code{e1071} needed), and
#' \code{'rf'} (random forest with default parameters, package \code{randomForest} needed).
#' If a list is used, it should contain two named components:
#' \code{fit} --- a function for model fitting, and
#' \code{predict} ---- a function for prediction.
#' 
#' @param train_percent The percentage (between 0 and 1) of labels that are used for model 
#' training during random splits
#' 
#' @param num_splits The number of random splits.
#' 
#' @param start_seed in the i-th split, the seed is set to start_seed + i.
#' 
#' @param verbose print progress every verbose splits if verbose is positive,
#' or remain quiet if verbose is zero
#' 
#' @return An object of class \code{PhecapModel}, with components
#' \item{model}{the fitted model}
#' \item{method}{the method used for model training}
#' \item{feature_selected}{the feature selected by SAFE}
#' \item{train_roc}{ROC on training dataset}
#' \item{train_auc}{AUC on training dataset}
#' \item{split_roc}{average ROC on random splits of training dataset}
#' \item{split_auc}{average AUC on random splits of training dataset}
#' \item{fit_function}{the function used for fitting}
#' \item{predict_function}{the function used for prediction}
#' 
#' @export
phecap_train_phenotyping_model <- function(
  data, surrogates, feature_selected, 
  method = "lasso_bic",
  train_percent = 0.7, num_splits = 200L,
  start_seed = 12345L, verbose = 10L)
{
  if (is.null(data$training_label)) {
    stop("Missing training label")
  }
  
  frame <- merge(
    data$icd_feature, data$nlp_feature, 
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$hu_feature,
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$training_label,
    by = data$patient_index, all = TRUE)

  feature <- phecap_generate_feature_matrix(
    data, surrogates, frame, feature_selected)
  label <- frame[, setdiff(
    names(data$training_label), data$patient_index)]
  ii <- which(!is.na(label))
  x <- feature[ii, , drop = FALSE]
  y <- label[ii]
  penalty_weight <- c(
    rep.int(0.0, attr(feature, "free")),
    rep.int(1.0, ncol(feature) - attr(feature, "free")))
  
  result <- get_roc_auc_with_splits(
    x, y, penalty_weight, method = method,
    train_percent = train_percent, num_splits = num_splits,
    start_seed = start_seed, verbose = verbose)
  if (is.numeric(result$model)) {
    names(result$model) <- c("(Intercept)", colnames(x))
  }
  result$feature_selected <- feature_selected
  
  class(result) <- "PhecapModel"
  return(result)
}


#' Validate the Phenotyping Model using the Validation Labels
#' 
#' Apply the trained model to all patients in the validation dataset,
#' and measure the prediction accuracy via ROC and AUC.
#' 
#' @usage phecap_validate_phenotyping_model(data, surrogates, model)
#' 
#' @param data an object of class \code{PhecapData}, 
#' obtained by calling \code{PhecapData(...)}
#' @param surrogates a list of objects of class \code{PhecapSurrogate}, obtained by 
#' something like \code{list(PhecapSurrogate(...), PhecapSurrogate(...))}.
#' @param model an object of class \code{PhecapModel}, obtained by calling 
#' \code{phecap_train_phenotyping_model}.
#' 
#' @return An object of class \code{PhecapValidation}, with components
#' \item{method}{the method used for model training}
#' \item{train_roc}{ROC on training dataset}
#' \item{train_auc}{AUC on training dataset}
#' \item{split_roc}{average ROC on random splits of training dataset}
#' \item{split_auc}{average AUC on random splits of training dataset}
#' \item{valid_roc}{ROC on validation dataset}
#' \item{valid_auc}{AUC on validation dataset}
#' 
#' @export
phecap_validate_phenotyping_model <- function(
  data, surrogates, model)
{
  if (is.null(data$validation_label)) {
    stop("Missing validation label")
  }
  
  frame <- merge(
    data$icd_feature, data$nlp_feature, 
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$hu_feature,
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$validation_label,
    by = data$patient_index, all = TRUE)

  feature <- phecap_generate_feature_matrix(
    data, surrogates, frame, model$feature_selected)
  label <- frame[, setdiff(
    names(data$validation_label), data$patient_index)]
  ii <- which(!is.na(label))
  x <- feature[ii, , drop = FALSE]
  y <- label[ii]

  prediction <- model$predict_function(model$model, x)
  valid_roc <- get_roc(y, prediction)
  valid_auc <- get_auc(y, prediction)
  
  result <- list(
    model = model$model,
    method = model$method,
    train_roc = model$train_roc,
    train_auc = model$train_auc,
    split_roc = model$split_roc,
    split_auc = model$split_auc,
    valid_roc = valid_roc,
    valid_auc = valid_auc)
  class(result) <- "PhecapValidation"
  
  return(result)
}


#' Predict Phenotype
#' 
#' Compute predicted probability of having the phenotype for each patient in the dataset.
#' 
#' @usage phecap_predict_phenotype(data, surrogates, model)
#' 
#' @param data an object of class \code{PhecapData}, obtained by calling \code{PhecapData(...)}.
#' @param surrogates a list of objects of class \code{PhecapSurrogate}, obtained by something like
#' \code{list(PhecapSurrogate(...), PhecapSurrogate(...))}.
#' @param model an object of class \code{PhecapModel}, probably returned from 
#' \code{phecap_train_phenotyping_model}.
#' 
#' @return A \code{data.frame} with two columns: 
#' \item{patient_index}{patient identifier,} 
#' \item{prediction}{predicted phenotype.}
#' 
#' @export
phecap_predict_phenotype <- function(
  data, surrogates, model)
{
  frame <- merge(
    data$icd_feature, data$nlp_feature, 
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$hu_feature,
    by = data$patient_index, all = TRUE)
  frame <- merge(
    frame, data$validation_label,
    by = data$patient_index, all = TRUE)

  feature <- phecap_generate_feature_matrix(
    data, surrogates, frame, model$feature_selected)
  prediction <- model$predict_function(model$model, feature)
  
  result <- data.frame(
    frame[[data$patient_index]],
    prediction)
  names(result) <- c(data$patient_index, "prediction")
  
  return(result)
}


print.PhecapFeatureExtraction <- function(x, ...)
{
  cat("Feature(s) selected by",
      "surrogate-assisted feature extraction (SAFE)\n")
  print(as.character(x), ...)
}


print.PhecapModel <- function(x, ...)
{
  cat("Phenotyping model:\n")
  print(x$model, ...)
  cat("AUC on training data:", 
      format(x$train_auc, digits = 3L), "\n")
  cat("Average AUC on random splits:", 
      format(x$split_auc, digits = 3L), "\n")
}


print.PhecapModel <- function(x, ...)
{
  cat("Phenotyping model:\n")
  print(x$model, ...)
  cat("AUC on training data:", 
      format(x$train_auc, digits = 3L), "\n")
  cat("Average AUC on random splits:", 
      format(x$split_auc, digits = 3L), "\n")
}


print.PhecapValidation <- function(x, ...)
{
  cat("AUC on validation data:", 
      format(x$valid_auc, digits = 3L), "\n")
  cat("AUC on training data:", 
      format(x$train_auc, digits = 3L), "\n")
  cat("Average AUC on random splits:", 
      format(x$split_auc, digits = 3L), "\n")
}


#' Plot ROC and Related Curves for Phenotyping Models
#' 
#' Plot ROC-like curves to illustrate phenotyping accuracy.
#' 
#' @usage phecap_plot_roc_curves(x, axis_x = "1 - spec", axis_y = "sen",
#' what = c("training", "random-splits", "validation"), ggplot = TRUE, ...)
#' 
#' @param x either a single object of class PhecapModel or PhecapValidation (returned from 
#' \code{phecap_train_phenotyping_model} or \code{phecap_validate_phenotyping_model}), 
#' or a named list of such objects
#' 
#' @param axis_x 
#' an expression that leads to the \code{x} coordinate.
#' Recognized quantities include:
#' \code{cut} (probability cutoff),
#' \code{pct} (percent of predicted cases),
#' \code{acc} (accuracy),
#' \code{tpr} (true positive rate),
#' \code{fpr} (false positive rate),
#' \code{tnr} (true negative rate),
#' \code{ppv} (positive predictive value),
#' \code{fdr} (false discovery rate),
#' \code{npv} (negative predictive value),
#' \code{sen} (sensitivity),
#' \code{spec} (specificity),
#' \code{prec} (precision),
#' \code{rec} (recall),
#' \code{f1} (F1 score).
#' @param axis_y an expression that leads to the \code{y} coordinate. 
#' Recognized quantities are the same as those in \code{axis_x}.
#' @param what The curves to be included in the figure.
#' @param ggplot if TRUE and ggplot2 is installed, ggplot will be used for the figure.
#' Otherwise, the base R graphics functions will be used.
#' @param \dots arguments to be ignored.
#' 
#' @export
phecap_plot_roc_curves <- function(
  x, axis_x = "1 - spec", axis_y = "sen", 
  what = c("training", "random-splits", "validation"),
  ggplot = TRUE, ...)  {
  object <- x
  if (is(x, "PhecapModel") || is(x, "PhecapValidation")) {
    object <- list(x)
    names(object) <- deparse(substitute(x))
  } else if (!is.list(x)) {
    stop("Not a PhecapModel / PhecapValidation object or a list of them")
  } else if (is.null(names(x))) {
    stop("List should be named")
  }
  
  df <- vector("list", length(object) * 3L)
  ii <- 1L
  for (kk in seq_along(object)) {
    oo <- object[[kk]]
    for (ww in what) {
      if (ww == "training" && "train_roc" %in% names(oo)) {
        df[[ii]] <- data.frame(
          cut = oo$train_roc$cut,
          value_x = eval(parse(text = axis_x), oo$train_roc),
          value_y = eval(parse(text = axis_y), oo$train_roc),
          kk = names(object)[kk], ww = ww)
      } else if (ww == "random-splits" && "split_roc" %in% names(oo)) {
        df[[ii]] <- data.frame(
          cut = oo$split_roc$cut,
          value_x = eval(parse(text = axis_x), oo$split_roc),
          value_y = eval(parse(text = axis_y), oo$split_roc),
          kk = names(object)[kk], ww = ww)
      } else if (ww == "validation" && "valid_roc" %in% names(oo)) {
        df[[ii]] <- data.frame(
          cut = oo$valid_roc$cut,
          value_x = eval(parse(text = axis_x), oo$valid_roc),
          value_y = eval(parse(text = axis_y), oo$valid_roc),
          kk = names(object)[kk], ww = ww)
      }
      ii <- ii + 1L
    }
  }

  df <- do.call("rbind", df)
  df <- aggregate(cbind(cut, value_y) ~ kk + ww + value_x, df, 
                  max)
  df <- df[order(df$kk, df$ww, df$cut, df$value_x, df$value_y), ]
  if (length(unique(df$kk)) > 1L) {
    if (length(unique(df$ww)) > 1L) {
      df$type <- paste(df$kk, df$ww, sep = ":")
    } else {
      df$type <- df$kk
    }
  } else if (length(unique(df$ww)) > 1L) {
    df$type <- df$ww
  }
  
  if (ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    pp <- ggplot2::ggplot(df)
    if ("type" %in% names(df)) {
      pp <- pp + ggplot2::geom_path(ggplot2::aes(
        x = df$value_x, y = df$value_y, color = type))
    } else {
      pp <- pp + ggplot2::geom_path(ggplot2::aes(
        x = df$value_x, y = df$value_y))
    }
    pp <- pp + ggplot2::xlab(axis_x) + ggplot2::ylab(axis_y)
    print(pp)
  } else {
    plot(NULL, NULL, 
         xlim = range(df$value_x), ylim = range(df$value_y),
         xlab = axis_x, ylab = axis_y)
    if ("type" %in% names(df)) {
      col <- 1L
      for (type in unique(df$type)) {
        lines(df[df$type == type, "value_x"],
              df[df$type == type, "value_y"], col = col)
        col <- col + 1L
      }
    } else {
      lines(df[, "value_x"], df[, "value_y"])
    }
    legend("bottomright", 
           legend = unique(df$type),
           lty = 1,
           col = seq_along(unique(df$type)))
  }
}
