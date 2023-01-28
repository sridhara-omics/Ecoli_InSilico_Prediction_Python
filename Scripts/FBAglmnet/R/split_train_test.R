#' Logical factor return to split data into train and test sets
#'
#' @param x - input matrix
#' @param ratio - split ratio into train and test
#'
#' @return
#' @export
#'
#' @examples
split_train_test <- function(x, ratio){
  ## set seed for reproducibility
  set.seed(123)

  ## dplyr function
  sample <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE, prob=c(1- ratio, ratio))

  return(sample)

}
