#' Plot an FPCA object. 
#'
#' \code{plot.FPCA} is currently implemented as just plotting the diagnostics plots, the same as CreateDiagnosticsPlot.
#' @param x An FPCA class object returned by FPCA()
#' @param ... passed into CreateDiagnosticsPlot
#' @export
#' @rdname CreateDiagnosticsPlot
plot.FPCA <- function(x, ...) {
  CreateDiagnosticsPlot(x, ...)
}