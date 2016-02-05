SetDerOptions <- function(derOptns = list()) {
  if (is.null(derOptns)) {
    derOptns <- list()
  }

  derOptns$p <- ifelse (is.null(derOptns$p), 0, derOptns$p)
  derOptns$method <- ifelse (is.null(derOptns$method), 'EIG',
                            derOptns$method)
  derOptns$GCV <- ifelse (is.null(derOptns$GCV), FALSE, TRUE)
  derOptns$kernelType <-  ifelse(is.null(derOptns$kernelType), 'epan',
                                 derOptns$kernelType)
  if (is.null(derOptns$h)) {
    derOptns$h <- derOptns$p * 0.1
  }

  derOptns
}
