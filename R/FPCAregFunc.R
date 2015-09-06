#' Functional Principal Component Analysis Regression with Functional dependent variable
#' 
#' Functional regression for dense or sparse functional data. 
#' 
#' @param expVarScal  A data.frame holding scalar explanatory variables; NAs will be omitted internally.
#' @param expVarFunc  A list holding the FPCA objects for each of the functional explanatory variables.
#' @param depVar  An FPCA object
#' @param regressionType A string defining the type of regression to perform ('dense' or 'sparse'); (default : automatically determined based on 'depVar')
#' @param bwScalar The value of bandwidth to be used for all scalar/functional cross-covariances (default: automatically determined using GCV)
#' @param bwFunct The values of bandwiths to be used for all function/function cross-covariances (default: automatically determined using GCV)
#' @param ...  Additional arguments 
#' 
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#' \cite{Senturk, D., Nguyen, D.V. "Varying Coefficient Models for Sparse Noise-contaminated Longitudinal Data", Statistica Sinica 21(4), (2011): 1831-1856. (Sparse data)}
#' @export

FPCAregFunc <- function(depVar,  expVarScal = NULL, expVarFunc = NULL, regressionType = NULL, bwScalar = NULL, bwFunct = NULL){
  
  if ( is.null(regressionType)){
    regressionType = depVar$optns$dataType
  }
  
  # Define number of predictors P and allocate the space for the beta(t)'s
  P = ncol(expVarScal) + length(expVarFunc)
  BetaFunctions = matrix( rep(0, length(depVar$workGrid) * P), nrow = P)

  # Centred and scale and numerical values / If it is a 2-D factor make it 0/1 
  Zvariables <- as.data.frame(sapply( expVarScal, function(x) 
                     if(is.numeric(x)){scale(x)}else if(is.factor(x)){ as.numeric(x)-1 } ))
  
  if( regressionType == 'Dense' ){
    fittedCurves = fitted(depVar)
    if (!is.null(expVarFunc)){
      fittedPredCurvesList <- lapply( expVarFunc, fitted)
    } else {
      fittedPredCurvesList <- NULL
    } 
    for (i in 1:length(depVar$workGrid)) {
      fitValuesAtTimeT = fittedCurves[,i];
      for (j in 1:length(expVarFunc)){
        eval(parse(text=paste("Zvariables$funcPredictor", j, " = fittedPredCurvesList[[j]][,i]", sep='')))
      }
      BetaFunctions[,i] =  coef( lm(fitValuesAtTimeT ~ ., data = Zvariables))[2: (1 + P)]
    }
  } else if( regressionType == 'Sparse'){
    # Get the Crosscovariance YZ
    CCYZ = matrix( rep(0, length(depVar$workGrid) *  ncol(expVarScal)), ncol =  ncol(expVarScal))
    for (i in 1:ncol(Zvariables)){
      CCYZ[,i] = CrCovYZ (Ymu=depVar$mu, Lt = depVar$inputData$t, Ly = depVar$inputData$y, Z= Zvariables[,i], 
                                    support= depVar$workGrid, bw = bwScalar)$smoothedCC
    }
    # If you have not functional predictors you are done!
    if (is.null(expVarFunc)){
      BetaFunctions = solve( cov( Zvariables), t(CCYZ))
    
    # Deal with functional predictors
    } else {

      Q =  length(expVarFunc)
      L =  length(expVarFunc[[1]]$workGrid)

      # Construct the matrix CCXX auto- and cross- covariances for the functional predictors
      CCXX = matrix( rep(0, ( L* Q)^2 ), nrow = L* Q ) 
      for(i in 1:( Q -0)){
        if( !all.equal(expVarFunc[[i]]$workGrid, depVar$workGrid) ){
          stop('You need the same workGrid for the functional predictors') 
        }
      }
      for( i in 1:Q){
        for( j in i:Q){
          if(i == j){
            CCXX[ ((i-1)*L) + (1:L), ((j-1)*L) + (1:L) ] = expVarFunc[[j]]$fittedCov
          } else {
            tempCCXX =  CrCovYX(Ly1 = expVarFunc[[i]]$inputData$y, Ly2 = expVarFunc[[j]]$inputData$y,
                               Lt1 = expVarFunc[[i]]$inputData$t, Lt2 = expVarFunc[[j]]$inputData$t,
                               Ymu1= expVarFunc[[i]]$mu,  Ymu2= expVarFunc[[j]]$mu, bw1 = bwFunct[1], bw2 = bwFunct[2]);
            CCXX[ ((i-1)*L) + (1:L), ((j-1)*L) + (1:L) ] = tempCCXX$smoothedCC;
            # And the transpose
            CCXX[ ((j-1)*L) + (1:L), ((i-1)*L) + (1:L) ] = t(tempCCXX$smoothedCC);
          }
        }  
      }
      
      # Construct the matrix CCYX for the cross-covariance between the func. predictors and the func. response variable
      CCYX =  matrix( rep(0, (L* L* Q) ), ncol = L )
      for (i in 1:Q){
        tempCCYX = CrCovYX(Ly1 = expVarFunc[[i]]$inputData$y, Ly2 = depVar$inputData$y,
                                Lt1 = expVarFunc[[i]]$inputData$t, Lt2 = depVar$inputData$t,
                                Ymu1= expVarFunc[[i]]$mu,  Ymu2= depVar$mu,  bw1 = bwFunct[1], bw2 = bwFunct[2]);
        CCYX[((i-1)*L) + (1:L),] = t(tempCCYX$smoothedCC);
      }

      # Construct the matrix CCXZ for the cross-covariance between the functional and the scalar predictors 
      CCXZ =  matrix( rep(0, (L * Q * ncol(expVarScal) ) ), ncol = L )
      for (i in 1:Q){
        for (j in 1:ncol(expVarScal)){
          CCXZ[j+ (i-1)* ncol(expVarScal), ] = 
                     CrCovYZ (Ymu = expVarFunc[[i]]$mu, Lt = expVarFunc[[i]]$inputData$t, Ly = expVarFunc[[i]]$inputData$y, 
                              Z = Zvariables[,j],support = expVarFunc[[i]]$workGrid, bw = bwScalar)$smoothedCC
        }
      }

      for( i in 1:length(depVar$workGrid)) {
        covYZX = matrix(rep(0,P),1)
        covZX = matrix(rep(0,P^2),P)
# browser()
        covYZX[1:(P-Q)] = CCYZ[i,];
        covYZX[(P-Q+1):P] =CCYX[ i+ (0:(Q-1))*L, i] ; 
       
        covZX[1:(P-Q),  1:(P-Q)] = cov(Zvariables);
        covZX[1:(P-Q),  (P-Q+1):P] = t( matrix( CCXZ[, i], Q, byrow=TRUE))
        covZX[(P-Q+1):P,  1:(P-Q)] =    matrix( CCXZ[, i], Q, byrow=TRUE)
        for( j0 in 1:Q){
          for( j1 in 1:Q){
            covZX[ P-Q+j0, P-Q+j1] = CCXX[i + (j0-1)*L , i + (j1-1)*L]
          }
        }     
        BetaFunctions[,i] = solve( covZX, t(covYZX))
      }
    }
  }
  FRegObj <- list(betaFunctions = BetaFunctions)
  return(FRegObj)
}
