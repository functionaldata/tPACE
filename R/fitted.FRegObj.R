#' Fitted functional sample from FPCA regression object
#' 
#' Combine the zero-meaned fitted values and the interpolated mean to get the final values.
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA().   
#' @param fpcaObjder A object of class FPCAder returned by the function FPCAder(. 
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' fittedY <- fitted(res)
#' @references
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export


fitted.FregObj <-  function (object, fpcaObj, scaleZ = FALSE, expVarScal, expVarFunc, ...) {
  
  fRegObj <- object;
  fpcaObj <- fpcaObj;
 
  if ( (is.null(expVarFunc) && !is.null(expVarScal)) ){

    # Centred and scale and numerical values / If it is a 2-D factor make it 0/1 
    if( scaleZ == FALSE) {
      Zvariables <- as.data.frame(sapply( expVarScal, function(x) 
                    if(is.numeric(x)){(x)}else if(is.factor(x)){ as.numeric(x)-1 } ))
    } else {
      Zvariables <- as.data.frame(sapply( expVarScal, function(x)  
                    if(is.numeric(x)){scale(x)}else if(is.factor(x)){ as.numeric(x)-1 } ))
    }

    ZMFV = as.matrix(Zvariables) %*% fRegObj$betaFunctions
    IM = approx(x= fpcaObj$obsGrid, y=fpcaObj$mu, fpcaObj$workGrid)$y 

    return( t(apply( ZMFV, 1, function(x) x + IM))) 
  } else if(  (!is.null(expVarFunc) && !is.null(expVarScal)) ) {
    
    
    # Centred and scale and numerical values / If it is a 2-D factor make it 0/1 
    if( scaleZ == FALSE) {
      Zvariables <- as.data.frame(sapply( expVarScal, function(x) 
                    if(is.numeric(x)){(x)}else if(is.factor(x)){ as.numeric(x)-1 } ))
    } else {
      Zvariables <- as.data.frame(sapply( expVarScal, function(x)  
                    if(is.numeric(x)){scale(x)}else if(is.factor(x)){ as.numeric(x)-1 } ))
    }

    ZMFV = as.matrix(Zvariables) %*% fRegObj$betaFunctions
    IM = approx(x= fpcaObj$obsGrid, y=fpcaObj$mu, fpcaObj$workGrid)$y 

    return( t(apply( ZMFV, 1, function(x) x + IM))) 
  }
}
