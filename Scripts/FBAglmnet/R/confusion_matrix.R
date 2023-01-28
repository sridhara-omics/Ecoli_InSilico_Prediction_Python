#' Confusion matrix - y versus yhat
#'
#' @param predicted_sources - predicted growth condition (encoded value)
#' @param actual_sources - actual growth condition (encoded value)
#' @param number_sources - total number of unique targets (including the train and test data)
#'
#' @return
#' @export
#'
#' @examples
confusion_matrix <- function(predicted_sources, actual_sources, number_sources){

  # fill out zeros for missing sources in the actual and predicted list
  pred_levels<-factor(as.numeric(paste(predicted_sources)),levels=c(1:number_sources))
  actual_levels<-factor(as.numeric(paste(actual_sources)),levels=c(1:number_sources))

  # confusion matrix to show actual values as columns and predicted as rows
  table_matrix<-table("Predicted"=pred_levels,"Actual"= actual_levels)

  return(table_matrix)

}
