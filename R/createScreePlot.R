#' Create the scree plot for the fitted eigenvalues
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param ys a vector of FVE (functional variation explained) from FPCA object
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createScreePlot(res$FVE)
#' @export

createScreePlot <-function(ys){ 
 
  dfbar <-barplot( rep(NA,length(ys)), ylim=c(0,105) ,axes=FALSE, xlab='Number of components', ylab='Fraction of Variance Explained', names.arg=1:length(ys), main="Screeplot")


  abline(h=(seq(0,100,5)), col="lightgray", lty="dotted")
  barplot(c(ys[1], diff(ys)), add = TRUE )
  lines(dfbar, y= ys, col='red')
  points(dfbar, y= ys, col='red')
  legend("right", "Cummul. FVE", col='red', lty=1, pch=1, bty='n') 
   
}
