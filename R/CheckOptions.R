#' Check option format
#'
#' Check if the options structure is valid and set the ones that are NULL
#' 
#' @param t is a n-by-1 list of vectors 
#' @param optns is an initialized option list
#' @param n is a total number of sample curves
#' @return logical

CheckOptions = function(t,optns,n){
  
    
  if( !(  any(optns[['useBinnedData']] == c('FORCE','AUTO','OFF')) )){ 
    # Force, turn off or automatically decide about the use of bin data
    cat("Error: FPCA is aborted because the argument: useBinnedData is invalid!\n"); 
    return(TRUE);   
  }
  if(  !( (length(optns[['bwmu']])==1) &&  is.numeric(optns[['bwmu']]) && (0<=optns[['bwmu']]) ) ){ 
    # bandwidth Bhoice for mean function is using CV or GCV
    cat("Error: FPCA is aborted because the argument: bwmu is invalid!\n"); 
    return(TRUE);   
  }
  if( !(  any(optns[['bwmuMethod']] == c('CV','GCV','GMeanAndGCV')) )){ 
    # bandwidth choice for mean function is GCV if bwmu = 0
    cat("Error: FPCA is aborted because the argument: bwmuMethod is invalid!\n"); 
    return(TRUE);   
  }
  if(!(length(optns[['bwuserCov']])==1) &&  is.numeric(optns[['bwuserCov']]) && (all(optns[['bwuserCov']]>=0))){ 
    # bandwidth choice for covariance function is CV or GCV
    cat("Error: FPCA is aborted because the argument: bwuserCov is invalid!\n"); 
    return(TRUE);   
  }
  if( !(  any(optns[['bwuserCovGcv']] == c('CV','GCV','GMeanAndGCV') ) )){ 
    # bandwidth choice for covariance function is GCV if bwuserCov = c(0,0)
    cat("Error: FPCA is aborted because the argument: bwuserCovGcv is invalid!\n");  
    return(TRUE);   
  }

  if (is.nan(optns[['kFoldCov']]) || optns[['kFoldCov']] < 2) {
    stop('Invalid `kFoldCov` option')
  }

  # if( !( (length(optns$ngrid1)==1) &&  is.numeric(optns$ngrid1) && (1<=optns$ngrid1) && (optns$ngrid1<90) ) ){ 
    # number of support points for the covariance surface 
  #   cat("Error: FPCA is aborted because the argument: ngrid1 is invalid!\n");  
  #   return(TRUE);   
  # }
  if( !(any(optns[['selectionMethod']] == c('FVE','AIC','BIC','fixedK')))){
    if ( !( is.numeric(optns[['selectionMethod']]) &&  (length(optns[['selectionMethod']])==1) && (1>=optns[['selectionMethod']]) && (optns[['selectionMethod']]<n) )){          
      # the method of choosing the number of principal components K
      cat("Error: FPCA is aborted because the argument: selectionMethod is invalid!\n");  
      return(TRUE);   
    }
  }
  if(  ( (length(optns[['FVEthreshold']])==1) &&  is.numeric(optns[['FVEthreshold']]) ) ){
    if (!( (0<=optns[['FVEthreshold']]) && (optns[['FVEthreshold']]<=1) ) ){  
      # the Fraction-of-Variance-Explained
      cat("Error: FPCA is aborted because the argument: FVEthreshold is invalid!\n"); 
      return(TRUE);   
    } 
  }
  if( !( (length(optns[['maxK']])==1) &&  is.numeric(optns[['maxK']]) && (1<=optns[['maxK']]) && (optns[['maxK']]<=n) )){  
    # maximum number of principal components to consider
    cat("Error: FPCA is aborted because the argument: maxK is invalid!\n");   
    return(TRUE);   
  } 
  #if( !is.null(optns$numComponents) ) {
  #  if( !( (length(optns$numComponents)==1) &&  is.numeric(optns$numComponents) && (1<=optns$numComponents) && (optns$numComponents<=n) )){  
  #    # maximum number of principal components to return
  #    cat("Error: FPCA is aborted because the argument: numComponents is invalid!\n");   
  #    return(TRUE);   
  #  }
  #} 
  if( !( is.null(optns[['dataType']]) || any(optns[['dataType']]==c("Sparse","DenseWithMV","Dense","p>>n")) )){ 
    #do we have regualr or sparse functional data
    cat("Error: FPCA is aborted because the argument: dataType is invalid!\n");     
    return(TRUE);     
  }   
  if( ( is.null(optns[['dataType']])  )){ 
 cat("Erroblaasdlf")
    optns[['dataType']] = IsRegular(t)
  }   
  if(!is.logical(optns[['error']])){ 
    # error assumption with measurement error 
    cat("Error: FPCA is aborted because the error option is invalid!\n");   
    return(TRUE);   
  }
  if( !( (length(optns[['nRegGrid']])==1) &&  is.numeric(optns[['nRegGrid']]) && (1<=optns[['nRegGrid']]) && (optns[['nRegGrid']]>optns[['maxK']]) ) ){
    # number of support points in each direction of covariance surface  
    cat("Error: FPCA is aborted because the argument: nRegGrid is invalid!\n");    
    return(TRUE);     
  }
  if( !(any(optns[['methodXi']] == c('CE','IN')))){ 
    #method to estimate the PC scores
    cat("Error: FPCA is aborted because the argument: methodXi is invalid!\n");   
    return(TRUE);   
  }
  #if(!is.logical(optns$shrink)){ 
  #  # apply shrinkage to estimates of random coefficients (dataType data only)
  #  cat("Error: FPCA is aborted because the argument: shrink is invalid!\n");   
  #  return(TRUE);   
  #}
  #if (! (  is.null(optns$newdata) || (is.numeric(optns$newdata) && is.vector(optns$newdata)))){
  #  # new time vector to evaluate the final estimates
  #  cat("Error: FPCA is aborted because the argument: newdata is invalid!\n");       
  #  return(TRUE);     
  #}
  if(!(any(optns[['kernel']] == c('epan','gauss','rect','quar','gausvar')))){ 
    #method to estimate the PC scores
    cat("Error: FPCA is aborted because the argument: kernel is invalid!\n");   
    return(TRUE);   
  }
  if( !( ( is.numeric(optns[['numBins']]) && (optns[['numBins']]>1)) || is.null(optns[['numBins']]) )  ){  
    # Check suitability of number of bins
    cat("Error: FPCA is aborted because the argument: numBins is invalid!\n");   
    return(TRUE);       
  }
  if( ( ( optns[['useBinnedData']] == 'FORCE') &&  is.null(optns[['numBins']]) ) ){  
    # Check that we have a number of the bins if we force binning
    cat("Error: FPCA is aborted because the argument: numBins is NULL but you FORCE binning!\n");   
    return(TRUE);       
  }
  if(!is.character(optns[['yname']])){ 
    # name of the variable analysed     
    cat("Error: FPCA is aborted because the argument: yname is invalid!\n");  
    return(TRUE);        
  }
  # if(!is.logical(optns$screePlot)){ 
  #   # make screeplot 
  #   cat("Error: FPCA is aborted because the argument: screePlot is invalid!\n");  
  #   return(TRUE);      
  # }
  # if(!is.logical(optns$designPlot)){ 
  #   # make designplot 
  #   cat("Error: FPCA is aborted because the argument: designPlot is invalid!\n");    
  #   return(TRUE);   
  # }
  # if(!is.logical(optns$corrPlot)){ 
  #   cat(optns$corrPlot)
  #   # make correlation plot 
  #   cat("Error: FPCA is aborted because the argument: corrPlot is invalid!\n");   
  #   return(TRUE);    
  # }
  if(!is.logical(optns[['diagnosticsPlot']])){ 
    # make diagnosticsPlot 
    cat("Error: FPCA is aborted because the argument: diagnosticsPlot is invalid!\n");    
    return(TRUE);   
  }
  if(!(any(optns[['rho']] == c('cv-random', 'cv', 'none', 'no')))){ 
    # truncation threshold for the iterative residual that is used 
    cat("Error: FPCA is aborted because the argument: rho is invalid!\n");     
    return(TRUE);   
  }
  if(!is.logical(optns[['verbose']])){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: verbose is invalid!\n");    
    return(TRUE);    
  }
  if (! (  is.null(optns[['userMu']]) || 
        (is.list(optns[['userMu']]) && is.vector(optns[['userMu']][['t']]) &&  is.vector(optns[['userMu']][['mu']]) &&
         ( length(optns[['userMu']][['t']]) ==  length(optns[['userMu']][['mu']]) ) ))){      
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: userMu is invalid!\n");     
    return(TRUE);   
  }

  if (! ( is.null(optns[['userCov']]) ||
        ( is.list(optns[['userCov']]) && is.vector(optns[['userCov']][['t']]) &&  is.matrix(optns[['userCov']][['cov']]) &&
          (length(optns[['userCov']][['t']]) ==  ncol(optns[['userCov']][['cov']]) ) && ( isSymmetric(optns[['userCov']][['cov']]) ) ) ) ){
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: userCov is invalid! (eg. Check if 'cov' is symmetric and 't' is of appropriate size.)\n");
    return(TRUE);
  }
  if(!(any(optns[['methodMu']] == c('PACE','RARE','CrossSectional')))){ 
    # user-defined mean functions
    cat("Error: FPCA is aborted because the argument: methodMu is invalid!\n");     
    return(TRUE);   
  }
  if( !( (length(optns[['outPercent']])==2) &&  is.numeric(optns[['outPercent']]) && all(0<=optns[['outPercent']]) && all(optns[['outPercent']]<=1) )){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: outPercent is invalid!\n");    
    return(TRUE);    
  }
  if( !( (length(optns[['rotationCut']])==2) &&  is.numeric(optns[['rotationCut']]) && all(0<=optns[['rotationCut']]) && all(optns[['rotationCut']]<=1) )){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: rotationCut is invalid!\n");    
    return(TRUE);    
  }
  if(is.logical(optns[['userCov']])){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: userCov is invalid!\n");     
    return(TRUE);   
  }
  
  
  return(FALSE);  
}

