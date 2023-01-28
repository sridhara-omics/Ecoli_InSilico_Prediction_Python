#' K-fold cross validation of glmnet
#'
#' @param inputs - features
#' @param targets - response variable
#' @param family - family of models e.g., binomial, multinomial, cox
#' @param alpha - elastic net, ridge regression or LASSO (lasso if set to 1.0)
#' @param nfolds - folds for cross validation
#'
#' @return
#' @export
#'
#' @examples
runGLMNET <- function(inputs, targets, family, alpha, nfolds) {
  inputs = as.matrix(inputs)
  targets = as.factor(targets)

  cvGLMnet<-cv.glmnet(inputs, targets, family=family, alpha=alpha, standardize=TRUE, maxit=900000, type.measure="class", nfolds=nfolds, lambda.min.ratio=0.0002, nlambda=2000)
  return(cvGLMnet)
}
