#' Create the scree plot for the fitted eigenvalues
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param fpcaObj A object of class FPCA returned by the function FPCA().
#' @param titleString a string variable to be used as title
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createScreePlot(res$cumFVE)
#' @export

createScreePlot <-function(fpcaObj, titleString = NULL){ 
 
  if(is.null(titleString)){
    titleString = "Screeplot"
  }
  ys <- fpcaObj$cumFVE;
  
  
  if( !is.vector(ys) ){ 
    stop('Please use a vector as input.')   
  }
  if(max(ys) > 100){
    warning('The maximum number in the input vector is larger than 100; are sure it is right?');
  }
  if(any(ys < 0) || any(diff(ys) <0) ){
    stop('This does not appear to be a valid cummulative FVE vector. Please check it carefully.')
  }

  dfbar <-barplot( rep(NA,length(ys)), ylim=c(0,105) ,axes=FALSE, 
                   xlab='Number of components', ylab='Fraction of Variance Explained', names.arg=1:length(ys), main=titleString)
  
  abline(h=(seq(0,100,5)), col="lightgray", lty="dotted")
  barplot(c(ys[1], diff(ys)), add = TRUE )
  lines(dfbar, y= ys, col='red')
  points(dfbar, y= ys, col='red')
  legend("right", "Cummul. FVE", col='red', lty=1, pch=1, bty='n') 
   
}
