#' Functional Principal Component Analysis Regression with Scalar dependent variable
#' 
#' Functional regression for dense or sparse functional data. 
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA(). 
#' @param extVar  A data.frame holding external explanatory variables.
#' @param depVar  A vector with the dependant variable.
#' @param varSelect  A string defining the type of step-wise variable selection applied ('AIC' or 'BIC'); this calls 'MASS::stepAIC()'. (default: NULL)
#' @param crossVal  A logical variable indicating if cross-validation measures for the predictive accuracy should be returned; this calls 'caret::train()'. (default: FALSE)
#' @param regressionType A string defining the type of regression to perform ('dense' or 'sparse'); (default : automatically determined based on 'fpcaObj'
#' @param ...  Additional arguments 
#' 
#' @details If both 'varSelect' and 'crossVal' are used then the cross-validation is applied to the model returned by the variable selection procedure.
#' 
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#' @export


FPCAregScalar <-  function (fpcaObj, extVar = NULL, depVar, varSelect = NULL, crossVal = FALSE, regressionType = NULL, ...) {
 
  if ( is.null(regressionType)){
    regressionType = fpcaObj$optns$dataType
  }
  
  if ( regressionType == 'Dense'){ 
    # Make dataset to regress on
    n = length(fpcaObj$lambda)
    Xi = data.frame(fpcaObj$xiEst)
    names(Xi) =  as.vector(mapply(  paste ,  rep('Xi',n), 1:n, sep=''))
    if ( is.null(extVar) ){ 
      theData = Xi
    } else {
      theData = data.frame( Xi, extVar);
    }
    
    # perform multiple linear regression
    lmObject <- lm( depVar ~ . , data = theData)
 
    # apply variable selection and/or cross-validation diagnostics if requested
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
    cvObject <- NULL
    if ( crossVal ){
      trainControl <- caret::trainControl(method="cv", number=10)
      cvObject <- caret::train( as.formula(lmObject), data=theData, trControl=trainControl, method="lm")  
    }    
    return( list(lmObject = lmObject, cvObject = cvObject) )

  } else if ( regressionType == 'Sparse'){
    print('Sparse regression is not yet implemented, contact Pantelis!')
    return(NULL)

  } else {
    stop('Unknown regression type requested.')
    return(NULL)
  } 
} 
