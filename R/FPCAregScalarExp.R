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
#' @param bwYZ A scalar for bandwidth used among along function-to-scalar cross-covariance estimations
#' @param bwYX A scalar for bandwidth used among along function-to-function cross-covariance estimations
#' @param ...  Additional arguments 
#' 
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. 
  #' "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#' \cite{Senturk, D., Nguyen, D.V. "Varying Coefficient Models for Sparse Noise-contaminated Longitudinal Data", 
#' Statistica Sinica 21(4), (2011): 1831-1856. (Sparse data)}
#' @export


FPCAregScalarExp <-  function (fpcaObjList, extVar = NULL, depVar, varSelect = NULL, bootStrap = FALSE, 
                               regressionType = NULL, y = NULL, t = NULL, lambda = 1e-9, bwYZ = NULL, bwYX = NULL, ...) {
  
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
    # Convenient reminder: 
    #     Y: scal. response, Z: scal. predictor, X: funct. predictor
    
    #stop("The sparse case is not ready yet.")
    #return(NULL)
    
    if (is.null(fpcaObjList)){
      stop('You have no functional covariates; just use regular regression functions.')
    }
    
    if (is.null(extVar)){
      nZ = 0;
    } else {
      nZ =  ncol(extVar)
    }
    
    nX = length(fpcaObjList)
    ngrid =  length(fpcaObjList[[1]]$workGrid)
    
    # Check that all functional objects have the same work-grid
    if (nX >1){
      for (i in 2:nX){
        if ( !all.equal(  fpcaObjList[[1]]$workGrid,  fpcaObjList[[i]]$workGrid    )){
          stop('All functional covariates need to have the same support. (still)')
        }
      } 
    }
    
    # Define basic arrays
    KovXYZ = array( 0, dim = c( ngrid, nX+nZ) )
    KovXXZZ = array( 0, dim = c(nX + nZ, nX + nZ, ngrid) )
    betaFuncs = array( 0, dim= c(ngrid, nZ + nX) )
    
    # Deal with Z-variables first if they exist
    if( nZ != 0 ){
      
      # Cross-covariance with X-variables
      for (j in 1:nX){ # X-variable counter
        for (i in 1:nZ){ # Z-variable counter
          
          y <- fpcaObjList[[j]]$inputData$y 
          t <- fpcaObjList[[j]]$inputData$t  
          KovXXZZ[nX+i,j,] = approx(x = fpcaObjList[[1]]$obsGrid, xout = fpcaObjList[[1]]$workGrid, y = 
                                      CrCovYZ(Z = extVar[,i], Ly = y, Lt=t, Ymu = fpcaObjList[[1]]$mu, 
                                              bw = bwYZ)$smoothedCC)$y
          KovXXZZ[j,nX+i,] = KovXXZZ[nX+i,j,];
        }
      }
      # Auto-covariance of Z's
      KovXXZZ[(nX+1):(nX+nZ), (nX+1):(nX+nZ),] = cov(extVar) 
      # Cross-covariance between Y and Z's
      KovXYZ[,(1+nX):(nX+nZ)] =  cov(depVar, extVar); 
    } 
    
    # Deal with X-variables next
    # Cross-covariance with other X-variables
    for (j in 1:nX){ # X-variable counter
      for (i in j:nX){ # Z-variable counter
        if (i == j){
          myDiag = diag(fpcaObjList[[j]]$fittedCov)
        } else {
          
          y1 <- fpcaObjList[[j]]$inputData$y 
          t1 <- fpcaObjList[[j]]$inputData$t  
          
          y2 <- fpcaObjList[[i]]$inputData$y 
          t2 <- fpcaObjList[[i]]$inputData$t   
          myDiag = diag(CrCovYX( Ly1 = y1, Lt1 = t1, Ly2 = y2, Lt2 = t2, fast = TRUE, bw1 = bwYX, bw2 = bwYX,
                                 Ymu1 = fpcaObjList[[i]]$mu, Ymu2 = fpcaObjList[[j]]$mu)$smoothedCC)
        }
        KovXXZZ[i,j,] = myDiag 
        KovXXZZ[j,i,] = myDiag 
      }
      bb13 = 123;
      KovXYZ[,j] = approx(x = fpcaObjList[[j]]$obsGrid, xout = fpcaObjList[[j]]$workGrid, y = 
                            CrCovYZ(Z = depVar, Ly = fpcaObjList[[j]]$inputData$y, 
                                    Lt= fpcaObjList[[j]]$inputData$t,
                                    Ymu = fpcaObjList[[j]]$mu, bw= bwYZ)$smoothedCC)$y
    }
    
    for (j in 1:length(fpcaObjList[[1]]$workGrid)){
      betaFuncs[j,] = solve ( a =  KovXXZZ[,,j], b= KovXYZ[j,] )
    }
    
    # K_x = 2;
    # phi = fpcaObjList[[1]]$phi
    # b = rep(0,K_x) 
    # for (k in 1:K_x){
    #  b[k] = crossprod( fpcaObjS$xiEst[,k], y = a)/crossprod( fpcaObjS$xiEst[,k] )
    # }
    # betaFuncs[,1] = crossprod( t(phi[,1:K_x]), b)
    # Notice that this works only for a single covariate
    
    return( list(lmObject = NULL, betaFunc = betaFuncs ))
    
  } else {
    stop('Unknown regression type requested.')
    return(NULL)
  } 
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

getBetas = function(data,indx ){  
  # indx is the random indexes for the bootstrap sample
  nsmall = dim(data)[2]; 
  return( as.numeric( qr.solve( a = data[indx,c(2:nsmall)], b = data[indx,1])) )  
}
