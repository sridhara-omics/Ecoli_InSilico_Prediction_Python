#' Predict growth sources using trained model and test data
#'
#' @param trained_model - cv.glmnet model
#' @param input_test - test features data
#'
#' @return
#' @export
#'
#' @examples
predictGLMNET <- function(trained_model, input_test){

  # derive coefficients from model
  s <- trained_model$lambda.min
  glmnet_coef <-as.vector(coef(trained_model, s=s))

  # predict using test data
  predicted_sources <- predict(trained_model, newx=as.matrix(input_test), type='class', s=s)

  return(predicted_sources)

}
