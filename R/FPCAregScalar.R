#' Functional Principal Component Analysis Regression with Scalar dependent variable
#' 
#' Functional regression for dense or sparse functional data. 
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA(). 
#' @param extVar  A data.frame holding external explanatory variables.
#' @param depVar  A vector with the dependant variable.
#' @param varSelect  A string defining the type of step-wise variable selection applied ('AIC' or 'BIC'); this calls 'MASS::stepAIC()'. (default: NULL)
#' @param regressionType A string defining the type of regression to perform ('dense' or 'sparse'); (default : automatically determined based on 'fpcaObj'
#' @param y A list of \emph{n} vectors containing the observed values for each individual.
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param lambda A scalar for the ridge correction in the sparse regression 
#' @param ...  Additional arguments 
#' 
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#' \cite{Senturk, D., Nguyen, D.V. "Varying Coefficient Models for Sparse Noise-contaminated Longitudinal Data", Statistica Sinica 21(4), (2011): 1831-1856. (Sparse data)}
#' @export


FPCAregScalar <-  function (fpcaObj, extVar = NULL, depVar, varSelect = NULL, bootStrap = FALSE, 
                            regressionType = NULL, y = NULL, t = NULL, lambda = 1e-9, ...) {
 
  if ( is.null(regressionType)){
    regressionType = fpcaObj$optns$dataType
  }
  
  if ( regressionType == 'Dense' ){ 
    
    # Use the FPCs to fit a functional linear regression using the 
    # decomposition of the functional regression into simple linear 
    # regressions between the FPC scores of the scalar response and
    # the predictors 
    
    # Impute the data in the case of sparse object and turn on varSelection
    if ( 'Sparse' == fpcaObj$optns$dataType ){ 
       imputedSample = fitted(fpcaObj);  
       s = fpcaObj$workGrid
       N = dim(fpcaObj$xiEst)[1];
       M = length(s);
       L3 = makePACEinputs(IDs = rep(1:N, each=M),tVec=rep(s,N), t(imputedSample) )
       fpcaObj = FPCA(y = L3$Ly, t = L3$Lt)
       if ( is.null(varSelect)){
         varSelect = 'AIC'
       }
    }

    # Make dataset to regress on
    n = length(fpcaObj$lambda)
    Xi = data.frame(fpcaObj$xiEst)
    names(Xi) =  as.vector(mapply( paste ,  rep('xiEstm',n), 1:n, sep=''))
    if ( is.null(extVar) ){ 
      theData = Xi
    } else {
      theData = data.frame( Xi, extVar);
    }
    
    # perform multiple linear regression
    lmObject <- lm( depVar ~ . , data = theData)
 
    # apply variable selection if requested
    if ( !is.null(varSelect)){
      if (varSelect == 'AIC') { 
        lmObject <- MASS::stepAIC(lmObject, trace = FALSE)     
      } else if ( varSelect == 'BIC') {
        lmObject <- MASS::stepAIC( lmObject, trace = FALSE, k = log(length(depVar)) )
      } else {  
        print("Invalid variable selection argument used; it must be 'AIC' or 'BIC'.")
       return(NULL)  
      }
    }
    
    # Get the Xi used
    coefNames <- names(coef(lmObject))
    xiEstIndcs = as.numeric(gsub("xiEstm", "",  coefNames[grep( "xiEstm",coefNames )]))      
    betaFunc = t(coef(lmObject)[grep( "xiEstm",coefNames )]) %*% t(fpcaObj$phi[ ,xiEstIndcs])
    
    bootBeta = NULL
    if (bootStrap){
      B = cbind( depVar, model.matrix(lmObject))
      bootBeta = boot::boot( data = B , statistic= getBetas, R = 1000 )      
    }
  
    return( list(lmObject = lmObject, betaFunc = as.vector(betaFunc), bootBeta =  bootBeta) )

  } else if ( regressionType == 'Sparse'){

    # Use the auto-covariance and cross-covariance matrices to fit a 
    # functional linear regression between the values of the scalar 
    # response and the predictors

    print('Sparse regression is not yet implemented, contact Pantelis!')
    return(NULL)

    # Check the data availability for further analysis
    if( is.null(y) || is.null(y) ){
      stop('Both y and t list have to be provided.')
      return(NULL)
    }    
    # Check the data validity for further analysis
    if( CheckData(y,t) ){
      cat('FPCA has stopped.')
      return(FALSE);
    }
    # Force the data to be list of numeric members
    y <- lapply(y, as.numeric)
    t <- lapply(t, as.numeric)

    # Convenient reminder: 
    #     Y: scal. response, Z: scal. predictor, X: funct. predictor
    
    if( !is.null(extVar) ){ 
      KKovYZ = cov(depVar, extVar); #Cross-covariance between y and z
      KovZZ = cov(extVar)
      KKovZX = CrCovYZ(Z = extVar, Ly = y, Lt=t)$smoothedCC # Semantic placeholder / this should fail
    } else { 
      KKovYZ = NULL
      KovZZ = NULL
      KKovZX = NULL
    }
           
    KovXX = fpcaObj$fittedCov
    KKovYX =  CrCovYZ(Z = depVar, Ly = y, Lt=t)$smoothedCC
    
    Knumer = makeCompositeCovMatrix(yx = KKovYX, yz = KKovYZ);
    Kdenom = makeCompositeCovMatrix(zz = KovZZ, xx = KovXX, zx = KovZX);
        
    betaFunc = solve(a = Kdenom + diag(x=lambda, nrow=nrow(Kdenom)), b = Knumer)
    return( lmObject = NULL, betaFunc = betaFunc )

  } else {
    stop('Unknown regression type requested.')
    return(NULL)
  } 
} 

makeCompositeCovMatrix <- function(yx = NULL, yz = NULL, xx = NULL, zz = NULL, zx = NULL){
  if(         is.null(yx) &&  is.null(yz) && !is.null(xx) &&  is.null(zz) &&  is.null(zx) ){
    K = xx
  } else if( !is.null(yx) &&  is.null(yz) &&  is.null(xx) &&  is.null(zz) &&  is.null(zx) ){
    K = yx
  } else if(  is.null(yx) &&  is.null(yz) && !is.null(xx) && !is.null(zz) && !is.null(zx) ){
    K =  matrix( c(zz, zx, rbind( zx, xx)  ),  dim(zz)[1] +  max(dim(zz)) )
  } else if(  is.null(yx) &&  is.null(yz) && !is.null(xx) && !is.null(zz) && !is.null(zx) ){
    K = rbind(yz, yx);
  } else {
    stop('Insufficient submatrices to construct composite matrix for regression.')    
    return(NULL)
  } 
  return(K)
}

getBetas = function(data,indx ){  
  # indx is the random indexes for the bootstrap sample
  nsmall = dim(data)[2]; 
  return( as.numeric( qr.solve( a = data[indx,c(2:nsmall)], b = data[indx,1])) )  
}

