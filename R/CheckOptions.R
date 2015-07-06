#' Check if the options structure is valid and set the ones that are NULL
#' 
#' @param n is a total number of sample curves
#' @param p is an initialized option list
#' @return logical
#' @examples 
#' 1 + 3

CheckOptions = function(t,optns,n){
  
  bwmu = optns$bwmu;                bwmuGcv = optns$bwmuGcv; 
  bwuserCov = optns$bwuserCov;            bwuserCovGcv = optns$bwuserCovGcv;
  ntest1 = optns$ntest1;            # ngrid1 = optns$ngrid1; 
  selectionMethod = optns$selectionMethod;  FVEthreshold = optns$FVEthreshold;
  maxK = optns$maxK;                
  dataType = optns$dataType;          error = optns$error; 
  nRegGrid = optns$nRegGrid;              method = optns$method; 
  shrink = optns$shrink;            newdata = optns$newdata; 
  kernel = optns$kernel;            numBins = optns$numBins; 
  yname = optns$yname;              screePlot = optns$screePlot; 
  designPlot = optns$designPlot;    rho = optns$rho;
  verbose = optns$verbose;          corrPlot = optns$corrPlot;
  userMu = optns$userMu;                  methodMu = optns$methodMu;
  outPercent = optns$outPercent;  userCov = optns$userCov
                                useBinnedData = optns$useBinnedData;
    
  if( !(  any(optns$useBinnedData == c('FORCE','AUTO','OFF')) )){ 
    # Force, turn off or automatically decide about the use of bin data
    cat("Error: FPCA is aborted because the argument: useBinnedData is invalid!\n"); 
    return(TRUE);   
  }
  if(  !( (length(optns$bwmu)==1) &&  is.numeric(optns$bwmu) && (0<=optns$bwmu) ) ){ 
    # bandwidth Bhoice for mean function is using CV or GCV
    cat("Error: FPCA is aborted because the argument: bwmu  is invalid!\n"); 
    return(TRUE);   
  }
  if( !(  any(optns$bwmuGcv == c('CV','GCV','GMeanAndGCV')) )){ 
    # bandwidth choice for mean function is GCV if bwmu = 0
    cat("Error: FPCA is aborted because the argument: bwmuGcv is invalid!\n"); 
    return(TRUE);   
  }
  if(!(length(optns$bwuserCov)==1) &&  is.numeric(optns$bwuserCov) && (all(optns$bwuserCov>=0))){ 
    # bandwidth choice for covariance function is CV or GCV
    cat("Error: FPCA is aborted because the argument: bwuserCov is invalid!\n"); 
    return(TRUE);   
  }
  if( !(  any(optns$bwuserCovGcv == c('CV','GCV','GMeanAndGCV') ) )){ 
    # bandwidth choice for covariance function is GCV if bwuserCov = c(0,0)
    cat("Error: FPCA is aborted because the argument: bwuserCovGcv is invalid!\n");  
    return(TRUE);   
  }
  if( !( (length(optns$ntest1)==1) &&  is.numeric(optns$ntest1) && (1<=optns$ntest1) && (optns$ntest1<n) )){ 
    # number of curves used for CV when choosing bandwidth  
    cat("Error: FPCA is aborted because the argument: ntest1 is invalid!\n");  
    return(TRUE);   
  }
  # if( !( (length(optns$ngrid1)==1) &&  is.numeric(optns$ngrid1) && (1<=optns$ngrid1) && (optns$ngrid1<90) ) ){ 
    # number of support points for the covariance surface 
  #   cat("Error: FPCA is aborted because the argument: ngrid1 is invalid!\n");  
  #   return(TRUE);   
  # }
  if( !(any(optns$selectionMethod == c('FVE','AIC','BIC')))){
    if ( !( is.numeric(optns$selectionMethod) &&  (length(optns$selectionMethod)==1) && (1>=optns$selectionMethod) && (optns$selectionMethod<n) )){          
      # the method of choosing the number of principal components K
      cat("Error: FPCA is aborted because the argument: selectionMethod is invalid!\n");  
      return(TRUE);   
    }
  }
  if(  ( (length(optns$FVEthreshold)==1) &&  is.numeric(optns$FVEthreshold) ) ){
    if (!( (0<=optns$FVEthreshold) && (optns$FVEthreshold<=1) ) ){  
      # the Fraction-of-Variance-Explained
      cat("Error: FPCA is aborted because the argument: FVEthreshold is invalid!\n"); 
      return(TRUE);   
    } 
  }
  if( !( (length(optns$maxK)==1) &&  is.numeric(optns$maxK) && (1<=optns$maxK) && (optns$maxK<=n) )){  
    # maximum number of principal components to consider
    cat("Error: FPCA is aborted because the argument: maxK is invalid!\n");   
    return(TRUE);   
  } 
  if( !( is.null(optns$dataType) || any(optns$dataType==c("Sparse","DenseWithMV","Dense","p>>n")) )){ 
    #do we have regualr or sparse functional data
    cat("Error: FPCA is aborted because the argument: dataType is invalid!\n");     
    return(TRUE);     
  }   
  if( ( is.null(optns$dataType)  )){ 
 cat("Erroblaasdlf")
    optns$dataType = IsRegular(t)
  }   
  if(!is.logical(optns$error)){ 
    # error assumption with measurement error 
    cat("Error: FPCA is aborted because the error option is invalid!\n");   
    return(TRUE);   
  }
  if( !( (length(optns$nRegGrid)==1) &&  is.numeric(optns$nRegGrid) && (1<=optns$nRegGrid) && (optns$nRegGrid>optns$maxK) ) ){
    # number of support points in each direction of covariance surface  
    cat("Error: FPCA is aborted because the argument: nRegGrid is invalid!\n");    
    return(TRUE);     
  }
  if( !(any(optns$method == c('CE','IN')))){ 
    #method to estimate the PC scores
    cat("Error: FPCA is aborted because the argument: method is invalid!\n");   
    return(TRUE);   
  }
  if(!is.logical(optns$shrink)){ 
    # apply shrinkage to estimates of random coefficients (dataType data only)
    cat("Error: FPCA is aborted because the argument: shrink is invalid!\n");   
    return(TRUE);   
  }
  if (! (  is.null(optns$newdata) || (is.numeric(optns$newdata) && is.vector(optns$newdata)))){
    # new time vector to evaluate the final estimates
    cat("Error: FPCA is aborted because the argument: newdata is invalid!\n");       
    return(TRUE);     
  }
  if(!(any(optns$kernel == c('epan','gauss','rect','quar','gausvar')))){ 
    #method to estimate the PC scores
    cat("Error: FPCA is aborted because the argument: kernel is invalid!\n");   
    return(TRUE);   
  }
  if( !( ( is.numeric(optns$numBins) && (optns$numBins>1)) || is.null(optns$numBins) )  ){  
    # Check suitability of number of bins
    cat("Error: FPCA is aborted because the argument: numBins is invalid!\n");   
    return(TRUE);       
  }
  if( ( ( optns$useBinnedData == 'FORCE') &&  is.null(optns$numBins) ) ){  
    # Check that we have a number of the bins if we force binning
    cat("Error: FPCA is aborted because the argument: numBins is NULL but you FORCE binning!\n");   
    return(TRUE);       
  }
  if(!is.character(optns$yname)){ 
    # name of the variable analysed     
    cat("Error: FPCA is aborted because the argument: yname is invalid!\n");  
    return(TRUE);        
  }
  if(!is.logical(optns$screePlot)){ 
    # make screeplot 
    cat("Error: FPCA is aborted because the argument: screePlot is invalid!\n");  
    return(TRUE);      
  }
  if(!is.logical(optns$designPlot)){ 
    # make designplot 
    cat("Error: FPCA is aborted because the argument: designPlot is invalid!\n");    
    return(TRUE);   
  }
  if(!is.logical(optns$corrPlot)){ 
    cat(optns$corrPlot)
    # make correlation plot 
    cat("Error: FPCA is aborted because the argument: corrPlot is invalid!\n");   
    return(TRUE);    
  }
  if(!(any(optns$rho == c('cv-random', 'cv', 'none', 'no')))){ 
    # truncation threshold for the iterative residual that is used 
    cat("Error: FPCA is aborted because the argument: rho is invalid!\n");     
    return(TRUE);   
  }
  if(!is.logical(optns$verbose)){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: verbose is invalid!\n");    
    return(TRUE);    
  }
  if (! (  is.null(optns$userMu) || (is.numeric(optns$userMu) && is.vector(optns$userMu)))){      
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: userMu is invalid!\n");     
    return(TRUE);   
  }
  if(!(any(optns$methodMu == c('PACE','RARE','CrossSectional')))){ 
    # user-defined mean functions
    cat("Error: FPCA is aborted because the argument: methodMu is invalid!\n");     
    return(TRUE);   
  }
  if( !( (length(optns$outPercent)==2) &&  is.numeric(optns$outPercent) && all(0<=optns$outPercent) && all(optns$outPercent<=1) )){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: outPercent is invalid!\n");    
    return(TRUE);    
  }
  if( !( (length(optns$rotationCut)==2) &&  is.numeric(optns$rotationCut) && all(0<=optns$rotationCut) && all(optns$rotationCut<=1) )){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: rotationCut is invalid!\n");    
    return(TRUE);    
  }
  if(is.logical(optns$userCov)){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: userCov is invalid!\n");     
    return(TRUE);   
  }
  
  
  return(FALSE);  
}

