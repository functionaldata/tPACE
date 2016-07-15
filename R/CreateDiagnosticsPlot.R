#' Functional Principal Component Analysis Diagnostics plot
#' 
#' Deprecated. Use \code{plot.FPCA} instead.
#' @param ... passed into \code{plot.FPCA}.
#' @export
#' @rdname plot.FPCA

CreateDiagnosticsPlot <-function(...){ 
  message('Deprecated and will be removed in the next version. Use plot.FPCA() instead.')
  plot.FPCA(...)
}

