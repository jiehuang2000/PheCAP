#' Phenotyping with EHR using CAP
#' 
#' Phenotyping with Eletronic Health Records using a Common Automated Pipeline
#' 
#' PheCAP provides a straightforward interface for conducting phenotyping on 
#' eletronic health records. One can specify the data via \code{\link{PhecapData}}, 
#' define surrogate using \code{\link{PhecapSurrogate}}. Next, one may run surrogate-assisted 
#' feature extraction (SAFE) by calling \code{\link{phecap_run_feature_extraction}}, and then 
#' train and validate phenotyping models via \code{\link{phecap_train_phenotyping_model}} and
#' \code{\link{phecap_validate_phenotyping_model}}. The predictive performance can be 
#' visualized using \code{\link{phecap_plot_roc_curves}}. Predicted phenotype is provided by 
#' \code{\link{phecap_predict_phenotype}}.
#'  
#' @docType package
#' 
#' @name PheCAP-package
#' 
#' @aliases PheCAP
NULL

#' @import glmnet
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom methods is
#' @importFrom stats aggregate
#' @importFrom stats approxfun
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats deviance
#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom stats predict
#' @importFrom stats quasibinomial
#' @importFrom utils read.csv
#' @importFrom utils read.table
NULL

