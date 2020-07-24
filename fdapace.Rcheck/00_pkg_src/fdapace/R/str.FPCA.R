#' Compactly display the structure of an FPCA object
#'
#' @param object An FPCA object
#' @param ... Not used
#'
#' @export
str.FPCA <- function(object, ...) {
  fpcaObj <- object
  NextMethod(max.level=1)
}
