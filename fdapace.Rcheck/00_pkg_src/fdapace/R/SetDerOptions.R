SetDerOptions <- function(fpcaObject = NULL, derOptns = list()) {
  if (is.null(derOptns)) {
    derOptns <- list()
  }
  # These are relevant for fitted.FPCA
  derOptns$method <- ifelse (is.null(derOptns$method), 'FPC',
                            derOptns$method)
  #derOptns$k <- ifelse (is.null(derOptns$k), length(fpcaObject$lambda), derOptns$k)
  # derOptns$GCV <- ifelse (is.null(derOptns$GCV), FALSE, TRUE)
  
  derOptns$p <- ifelse (is.null(derOptns$p), 1, derOptns$p)
  derOptns$kernelType <-  ifelse(is.null(derOptns$kernelType), 'gauss',
                                 derOptns$kernelType)
  if (is.null(derOptns$bwMu) && is.null(derOptns$bwCov)) {
    if (is.null(derOptns$bw)) {
      derOptns$bw <- 
        if (!is.null(fpcaObject[['sigma2']]) && (fpcaObject$sigma2 / sum(fpcaObject$lambda)) >= 0.01) {
          derOptns$p * 0.10 * diff(range(fpcaObject$workGrid))
        } else {
          derOptns$p * 0.05 * diff(range(fpcaObject$workGrid))
        }
    }
    derOptns$bwCov <- derOptns$bwMu <- derOptns$bw
  } else if (!is.null(derOptns$bwMu) && !is.null(derOptns$bwCov)) {
    # OK
  } else {
    stop('need to specify neither or both bwMu and bwCov')
  }

  return(derOptns)
}
