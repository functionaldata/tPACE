#' Create estimates for the beta functions used
#'
#' Create estimates for the beta functions used in a FPCA regression object
#'
#' @param fpcaRegObj A object of class FPCAreg returned by the function FPCAregScalar().
#' @param sList A list with the support of each function / If vector then the same support is used for all betas.
#' 
#' @export


createBetaPlots <- function(fpcaRegObj, sList, openNewDev = TRUE){
 
  if(openNewDev){ 
    dev.new(width=7.95, height=3.0, noRStudioGD=TRUE) ; 
  }
  

  noOfBetas <- length(fpcaRegObj$betaFunc)  
  xlabelString = 's';  
  par(mfrow = c(1,noOfBetas))
  yedge = max(abs(unlist(fpcaRegObj$betaFunc))) *1.1
  
  if( is.vector(sList) &&  !is.list(sList) ){
    sList =rep(list(sList),noOfBetas)
  }
  
  for (i in 1:noOfBetas){
    ylabelString = paste('beta function ',  i ,   sep=''   )
    plot(sList[[i]], fpcaRegObj$betaFunc[[i]], type= 'l', 
      xlab= xlabelString, ylab=ylabelString    )
    abline(v=(seq(min(sList[[i]]), max(sList[[i]]), length.out=  min(30, diff( range(sList[[i]])) +1 ))), 
      col="lightgray", lty="dotted")
  }
}



