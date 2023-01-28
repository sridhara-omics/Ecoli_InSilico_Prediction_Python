#' Misclassification rate - yhat not predicted as y
#'
#' @param confusion_table - actual versus predicted as a table
#'
#' @return
#' @export
#'
#' @examples
misclassification_rate <- function(confusion_table) {

  rate <- 100 * (sum(confusion_table)-sum(diag(confusion_table))) / sum(confusion_table)

  return(rate)

}
