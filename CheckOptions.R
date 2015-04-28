CheckOptions = function(p,n){
  
  # Check if the options structure is valid
  # n : total number of sample curves
  # p : initialized option list
  
  bwmu = p$bwmu;                bwmu_gcv = p$bwmu_gcv; 
  bwxcov = p$bwxcov;            bwxcov_gcv = p$bwxcov_gcv;
  ntest1 = p$ntest1;            ngrid1 = p$ngrid1; 
  selection_k = p$selection_k;  FVE_threshold = p$FVE_threshold;
  maxk = p$maxk;                
  regular = p$regular;          error = p$error; 
  ngrid = p$ngrid;              method = p$method; 
  shrink = p$shrink;            newdata = p$newdata; 
  kernel = p$kernel;            numBins = p$numBins; 
  yname = p$yname;              screePlot = p$screePlot; 
  designPlot = p$designPlot;    rho = p$rho;
  verbose = p$verbose;          corrPlot = p$corrPlot;
  xmu = p$xmu;                  method_mu = p$method_mu;
  out_percent = p$out_percent;  xcov = p$xcov
  
  
  if(  !( (length(p$bwmu)==1) &  is.numeric(p$bwmu) & (0>=p$bwmu) ) ){ 
    # bandwidth choice for mean function is using CV or GCV
    cat("Error: FPCA is aborted because the argument: bwmu  is invalid!\n"); 
    return(TRUE);   
  }
  if( !(  any(p$bwmu_gcv == c(1,2,3)) )){ 
    # bandwidth choice for mean function is GCV if bwmu = 0
    cat("Error: FPCA is aborted because the argument: bwmu_gcv is invalid!\n"); 
    return(TRUE);   
  }
  if(!(length(p$bwxcov)==2) &  is.numeric(p$bwxcov) & (all(p$bwxcov>=0))){ 
    # bandwidth choice for covariance function is CV or GCV
    cat("Error: FPCA is aborted because the argument: bwxcov is invalid!\n"); 
    return(TRUE);   
  }
  if( !(  any(p$bwxcov_gcv == c(1,2,3)) )){ 
    # bandwidth choice for covariance function is GCV if bwxcov = c(0,0)
    cat("Error: FPCA is aborted because the argument: bwxcov_gcv is invalid!\n");  
    return(TRUE);   
  }
  if( !( (length(p$ntest1)==1) &  is.numeric(p$ntest1) & (1<=p$ntest1) & (p$ntest1<n) )){ 
    # number of curves used for CV when choosing bandwidth  
    cat("Error: FPCA is aborted because the argument: ntest1 is invalid!\n");  
    return(TRUE);   
  }
  if( !( (length(p$ngrid1)==1) &  is.numeric(p$ngrid1) & (1<=p$ngrid1) & (p$ngrid1<90) ) ){ 
    # number of support points for the covariance surface 
    cat("Error: FPCA is aborted because the argument: ngrid1 is invalid!\n");  
    return(TRUE);   
  }
  if( !(any(p$selection_k == c('FVE','AIC','BIC')))){
    if ( !( is.numeric(p$selection_k) &  (length(p$selection_k)==1) & (1>=p$selection_k) & (p$selection_k<n) )){          
      # the method of choosing the number of principal components K
      cat("Error: FPCA is aborted because the argument: selection_k is invalid!\n");  
      return(TRUE);   
    }
  }
  if( !( (length(p$FVE_threshold)==1) &  is.numeric(p$FVE_threshold) & (0<=p$FVE_threshold) & (p$FVE_threshold<=1) )){  
    # the Fraction-of-Variance-Explained
    cat("Error: FPCA is aborted because the argument: FVE_threshold is invalid!\n"); 
    return(TRUE);    
  }
  if( !( (length(p$maxk)==1) &  is.numeric(p$maxk) & (1<=p$maxk) & (p$maxk<=n) )){  
    # maximum number of principal components to consider
    cat("Error: FPCA is aborted because the argument: maxk is invalid!\n");   
    return(TRUE);   
  } 
  if( !( is.null(p$regular) || any(p$regular==c(1,2,3)) )){ 
    #do we have regualr or sparse functional data
    cat("Error: FPCA is aborted because the argument: regular is invalid!\n");     
    return(TRUE);     
  }   
  if(!is.logical(p$error)){ 
    # error assumption with measurement error 
    cat("Error: FPCA is aborted because the error option is invalid!\n");   
    return(TRUE);   
  }
  if( !( (length(p$ngrid)==1) &  is.numeric(p$ngrid) & (1<=p$ngrid) & (p$ngrid>p$maxk) ) ){
    # number of support points in each direction of covariance surface  
    cat("Error: FPCA is aborted because the argument: ngrid is invalid!\n");    
    return(TRUE);     
  }
  if( !(any(p$method == c('CE','IN')))){ 
    #method to estimate the PC scores
    cat("Error: FPCA is aborted because the argument: method is invalid!\n");   
    return(TRUE);   
  }
  if(!is.logical(p$shrink)){ 
    # apply shrinkage to estimates of random coefficients (regular data only)
    cat("Error: FPCA is aborted because the argument: shrink is invalid!\n");   
    return(TRUE);   
  }
  if (! (  is.null(p$newdata) || (is.numeric(p$newdata) & is.vector(p$newdata)))){
    # new time vector to evaluate the final estimates
    cat("Error: FPCA is aborted because the argument: newdata is invalid!\n");       
    return(TRUE);     
  }
  if(!(any(p$kernel == c('epan','gauss','rect','quar','gausvar')))){ 
    #method to estimate the PC scores
    cat("Error: FPCA is aborted because the argument: kernel is invalid!\n");   
    return(TRUE);   
  }
  if( !( is.numeric(p$numBins) || is.null(p$numBins) )  ){  
    # Check suitability of number of bins
    cat("Error: FPCA is aborted because the argument: numBins is invalid!\n");   
    return(TRUE);       
  }
  if(!is.character(p$yname)){ 
    # name of the variable analysed     
    cat("Error: FPCA is aborted because the argument: yname is invalid!\n");  
    return(TRUE);        
  }
  if(!is.logical(p$screePlot)){ 
    # make screeplot 
    cat("Error: FPCA is aborted because the argument: screePlot is invalid!\n");  
    return(TRUE);      
  }
  if(!is.logical(p$designPlot)){ 
    # make designplot 
    cat("Error: FPCA is aborted because the argument: designPlot is invalid!\n");    
    return(TRUE);   
  }
  if(!is.logical(p$corrPlot)){ 
    cat(p$corrPlot)
    # make correlation plot 
    cat("Error: FPCA is aborted because the argument: corrPlot is invalid!\n");   
    return(TRUE);    
  }
  if(!(any(p$rho == c('cv-random','cv','none')))){ 
    # truncation threshold for the iterative residual that is used 
    cat("Error: FPCA is aborted because the argument: rho is invalid!\n");     
    return(TRUE);   
  }
  if(!is.logical(p$verbose)){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: verbose is invalid!\n");    
    return(TRUE);    
  }
  if (! (  is.null(p$xmu) || (is.numeric(p$xmu) & is.vector(p$xmu)))){      
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: xmu is invalid!\n");     
    return(TRUE);   
  }
  if(!(any(p$method_mu == c('PACE','RARE')))){ 
    # user-defined mean functions
    cat("Error: FPCA is aborted because the argument: method_mu is invalid!\n");     
    return(TRUE);   
  }
  if( !( (length(p$out_percent)==1) &  is.numeric(p$out_percent) & (0<=p$out_percent) & (p$out_percent<=1) )){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: out_percent is invalid!\n");    
    return(TRUE);    
  }
  if(is.logical(p$xcov)){ 
    # display diagnostic messages
    cat("Error: FPCA is aborted because the argument: xcov is invalid!\n");     
    return(TRUE);   
  }
  
  
  return(FALSE);  
}

