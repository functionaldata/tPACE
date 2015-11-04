#' Functional Principal Component Analysis Regression with Scalar dependent variable
#' 
#' Functional regression for dense or sparse functional data. 
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA(). 
#' @param extVar  A data.frame holding external explanatory variables; NAs will be omitted internally.
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


FPCAregScalar <-  function (fpcaObjList, extVar = NULL, depVar, varSelect = NULL, bootStrap = FALSE, 
                            regressionType = NULL, y = NULL, t = NULL, lambda = 1e-9, ...) {
  
  #If it is a single element automatically coerce it to be a single element list
  if (class(fpcaObjList) == "FPCA"){
    fpcaObjList = list(fpcaObjList)
  }
  # Check that the list has only FPCA-class elements
  if( !all(sapply( fpcaObjList, class) == "FPCA") ){
    stop("Invalid fpcaObjList; all elements of it must be of class FPCA.")
    return(NULL)
  }
  # Check that the samples are comparable
  if( (length(fpcaObjList)  > 1) && ( diff(range(sapply( fpcaObjList, function(x) nrow(x[['xiEst']]))) ) != 0) ){
    stop("Invalid fpcaObjList; the xiEst have different length.")
    return(NULL)
  }
  
  if ( is.null(regressionType)){
    regressionType = fpcaObjList[[1]]$optns$dataType
  }
  
  if ( regressionType == 'Dense' ){ 
    
    # Use the FPCs to fit a functional linear regression using the 
    # decomposition of the functional regression into simple linear 
    # regressions between the FPC scores of the scalar response and
    # the functional predictors 
    
    # Impute the data in the case of sparse object and turn on varSelection
    for (j in 1:length(fpcaObjList) ){
      if ( 'Sparse' == fpcaObjList[[j]]$optns$dataType ){
        fpcaObjList[[j]] = makeDenseObj( fpcaObjList[[j]] );
        if ( is.null(varSelect)){
          varSelect = 'AIC'
        }
      }
    }
    
    # Make dataset to regress on
    plengths = sapply( fpcaObjList, function(x) ncol(x[['xiEst']]))
    p = sum( plengths );
    q = nrow( fpcaObjList[[1]][['xiEst']] )
    Xi =  data.frame(matrix( rep(0,q*p), ncol = p))
    # We set names equal to xiEst1A, xiEst2A,.. xiEstnA, xiEst1B, xiEst2B,....
    names(Xi) = unlist(mapply( function(ll,i) as.vector( mapply( paste,  
                                                                 rep('xiEstm',ll), LETTERS[i] ,1:ll, sep='')), plengths, seq_along(plengths) ))
    
    smallpEnd = cumsum(sapply( fpcaObjList, function(x) length(x[['lambda']])))   
    smallpBeg = c(1, smallpEnd + 1); smallpBeg = smallpBeg[-length(smallpBeg)] 
    for ( j in 1:length(fpcaObjList) ){
      Xi[ ,c(smallpBeg[j]:smallpEnd[j])] = fpcaObjList[[j]][['xiEst']]
    }
    if ( is.null(extVar) ){ 
      theData = data.frame( Xi, depVar);
    } else {
      theData = data.frame( Xi, extVar,depVar);
    }
    
    
    # perform multiple linear regression
    lmObject <- lm( depVar ~ . , data = na.omit(theData))
    
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
    betaFunc = list();
    
    for(j in 1:length(fpcaObjList)){
      xiNames <- paste("xiEstm",LETTERS[j],sep='')
      xiEstIndcs = as.numeric(gsub(xiNames, "",  coefNames[grep( xiNames,coefNames )]))  
      betaFunc[[j]] = as.vector( t(coef(lmObject)[grep( xiNames,coefNames )]) %*% t(fpcaObjList[[j]]$phi[ ,xiEstIndcs]) )
    } 
    
    bootBeta = NULL
    if (bootStrap){
      B = cbind( depVar, model.matrix(lmObject))
      bootBeta = boot::boot( data = B , statistic= getBetas, R = 1000 )      
    }
    
    return( list(lmObject = lmObject, betaFunc = betaFunc, bootBeta =  bootBeta) )
    
  } else if ( regressionType == 'Sparse'){
    
    # Use the auto-covariance and cross-covariance matrices to fit a 
    # functional linear regression between the values of the scalar 
    # response and the predictors
    
    print('Sparse regression is not yet implemented, contact Pantelis!')
    return(NULL)
    
    y <- fpcaObjList[[1]]$inputData$y 
    t <- fpcaObjList[[1]]$inputData$t
    
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


makeDenseObj = function( fo1 ){
  # Impute the data in the case of sparse object and turn on varSelection
  imputedSample = fitted(fo1);
  s = fo1$workGrid
  N = dim(fo1$xiEst)[1];
  M = length(s);
  L3 = makePACEinputs(tVec= s, yVec = imputedSample )
  fo2 = FPCA(y = L3$Ly, t = L3$Lt)
  return( fo2 )
}


