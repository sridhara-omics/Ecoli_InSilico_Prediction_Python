#' Set very small numbers to zero
#'
#' @param x - float matrix
#' @param eps - matrix values below threshold are set to 0
#'
#' @return
#' @export
#'
#' @examples
threshold <- function(x, eps) {
  x[abs(x) < eps]=0
  return(x)
}
