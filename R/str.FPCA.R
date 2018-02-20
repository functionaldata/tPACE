#' Compactly display the structure of an FPCA object
#'
#' object An FPCA object
#' ... Not used
#'
#' @export
str.FPCA <- function(object, ...) {
  fpcaObj <- object
  NextMethod(max.level=1)
}
