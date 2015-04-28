SetOptions = function(bwmu = 0, bwmu_gcv = 1, bwxcov = c(0,0), bwxcov_gcv = 1, 
    ntest1 = 30, ngrid1 = 30, selection_k = "BIC", FVE_threshold = 0.95,
    maxk = 20, regular = NULL, error = TRUE, ngrid = 51,
    method = "CE", shrink = FALSE, newdata = NULL, kernel = "gauss", 
    numBins = NULL, yname = NULL, screePlot = FALSE, designPlot = FALSE, 
    corrPlot = FALSE,     rho = "cv", verbose = TRUE, xmu = NULL, xcov = NULL, 
    method_mu = 'PACE', out_percent = 0){
  
  # Set the PCA option structure.
  # bwmu : bandwidth choice for mean function is using CV or GCV
  # bwmu_gcv : bandwidth choice for mean function is GCV if bwmu = 0
  # bwcov : bandwidth choice for covariance function is CV or GCV
  # bwcov_gcv : bandwidth choice for covariance function is GCV if bwxcov = c(0,0)
  # ntest1 : number of curves used for CV when choosing bandwidth 
  # ngrid1 : number of support points for the covariance surface 
  # selection_k : the method of choosing the number of principal components K
  # FVE_threshold : Fraction-of-Variance-Explained
  # maxk : maximum number of principal components to consider
  # regular : do we have regular or sparse functional data
  # error : error assumption with measurement error
  # ngrid : number of support points in each direction of covariance surface 
  # method : method to estimate the PC scores
  # shrink : apply shrinkage to estimates of random coefficients (regular data only)
  # newdata : new data points to estimate
  # kernel : smoothing kernel choice
  # numBins : number of bins
  # yname : name of the variable analysed
  # screePlot : make scree plot
  # designPlot : make design plot
  # corrPlot : make correlation plot
  # rho : truncation threshold for the iterative residual  
  # verbose : display diagnostic messages
  # xmu : user-defined smoothed mean function 
  # xcov : user-defined smoothed covariance function
  # method_mu :  method to estimate mu
  # out_percent : number in [0,1] indicating the out_percent data in the boundary
  
  
  if(is.null(bwmu)){ # bandwidth choice for mean function is using CV or GCV
    bwmu = 0;   
  }
  if(is.null(bwmu_gcv)){ # bandwidth choice for mean function is GCV if bwmu = 0
    bwmu_gcv = 1;  
  }
  if(is.null(bwxcov)){ # bandwidth choice for covariance function is CV or GCV
    bwxcov = c(0,0); 
  }
  if(is.null(bwxcov_gcv)){  # bandwidth choice for covariance function is GCV if bwxcov = c(0,0)
    bwxcov_gcv = 1;
  }
  if(is.null(ntest1)){ # number of curves used for CV when choosing bandwidth 
    ntest1 = 30;
  }
  if(is.null(ngrid1)){ # number of support points for the covariance surface 
    ngrid1 = 30;
  }
  if(is.null(selection_k)){ # the method of choosing the number of principal components K
    selection_k = "BIC";
  }
  if(selection_k == "FVE" && is.null(FVE_threshold)){  # the Fraction-of-Variance-Explained
    FVE_threshold = 0.95;
  }
  if(is.null(maxk)){ # maximum number of principal components to consider
    maxk = 20;   
  }
  if(is.null(regular)){ #do we have regular or sparse functional data
    regular = NULL;    
  }
  if(is.null(error)){ # error assumption with measurement error
    error = TRUE;    
  }
  if(is.null(ngrid)){ # number of support points in each direction of covariance surface 
    ngrid = 51;    
  }
  if(is.null(method)){ # method to estimate the PC scores
    shrink = FALSE;
    method = "CE";   
  }
  if(is.null(shrink)){ # apply shrinkage to estimates of random coefficients (regular data only)
    shrink = FALSE;
  }
  if(shrink == TRUE && (error != TRUE || method != "IN")){ # Check for valid shrinkage choice
    cat('shrinkage method only had effects when method = "IN" and error = TRUE! Reset to shrink = FALSE now!\n');
    shrink = FALSE      
  }
  if(is.null(newdata)){ # new data points to estimnate
    newdata <- NULL
  }
  if(is.null(kernel)){ # smoothing kernel choice
    if(regular == 2){
      kernel = "epan";   # kernel: Epanechnikov
    }else{
      kernel = "gauss";  # kernel: Gaussian
    }
  }
  if(is.null(yname)){ # name of the variable analysed
    yname = as.character(substitute(y))      
  }
  if(maxk > (ngrid-2)){ # check if a reasonable number of eigenfunctions is requested
    cat(paste("maxk can only be less than or equal to", ngrid-2,"! Reset to be", ngrid-2, "now!\n"));
    maxk = ngrid -2;
  }
  if(is.numeric(selection_k)){
    if(selection_k > (ngrid-2)){ # check if a reasonable number of eigenfunctions is requested
      cat(paste("maxk can only be less than or equal to", ngrid-2,"! Reset to be", ngrid-2, "now!\n"));
      maxk = ngrid -2;
    }else if(selection_k <= 0){ # check if a positive number of eigenfunctions is requested
      cat("selection_k must be a positive integer! Reset to BIC now!\n");
      selection_k = "BIC"
      FVE_threshold = 0.95;
    }
  }
  if(is.null(screePlot)){ # make screeplot
    screePlot = FALSE;
  }
  if(is.null(designPlot)){ # make designplot
    designPlot = FALSE;
  }
  if(is.null(corrPlot)){ # make corrplot
    corrPlot = FALSE;
  }
  if(is.null(rho)){ # truncation threshold for the iterative residual that is used
    rho = "cv";
  }
  if(is.null(verbose)){ # display diagnostic messages
    verbose = TRUE;
  }  
  if(is.null(newdata)){ # new data points to estimate
    newdata <- NULL
  }
  if(is.null(xmu)){ # user-defined mean functions valued at distinct input time points
    xmu <- NULL
  }
  if(is.null(xcov)){
    xcov <- NULL
  }
  if(is.null(method_mu)){ # method to estimate mu
    method_mu <- 'PACE'
  }
  if(is.null(out_percent)){ # number in [0,1] indicating if we leave out out_percent data in the boundary
    out_percent <- 0
  }  
  if(error == FALSE && (selection_k == "AIC" || selection_k == "BIC")){ # Check suitability of information criterion
    cat('When assume no measurement error, cannot use "AIC" or "BIC". Reset to "BIC" now!\n')
    selection_k = "BIC" 
  }
  kernNames = c("rect", "gauss", "epan", "gausvar", "quar");
  if(!(kernel %in% kernNames)){ # Check suitability of kernel
    cat(paste('kernel', kernel, 'is unrecognizable! Reset to Epanechnikov kernel now!\n')); 
    kernel = "epan"; 
  }  
  if(!is.null(numBins)){ 
    if(numBins < 10){   # Check suitability of number of bins
      cat("Number of bins must be at least +10! No binning will be performed!\n");
      numBins = 0;
    }
  }  
  return( list(bwmu = bwmu, bwmu_gcv = bwmu_gcv, bwxcov = bwxcov, bwxcov_gcv = bwxcov_gcv,
          ntest1 = ntest1, ngrid1 = ngrid1, selection_k = selection_k, FVE_threshold = FVE_threshold,
          maxk = maxk, regular = regular, error = error, ngrid = ngrid, 
          method = method, shrink = shrink, newdata = newdata, kernel = kernel, corrPlot = corrPlot,	
          numBins = numBins, yname = yname, screePlot = screePlot, designPlot = designPlot, rho = rho, 
          verbose = verbose, xmu = xmu, xcov = xcov, method_mu= method_mu, out_percent = out_percent) )
}
