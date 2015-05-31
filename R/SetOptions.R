#' Set the PCA option list
#' 
#' @param bwmu : bandwidth choice for mean function is using CV or GCV
#' @param bwmu_gcv : bandwidth choice method for mean function is GCV if bwmu = 0
#' @param bwcov : bandwidth choice for covariance function is CV or GCV
#' @param bwcov_gcv : bandwidth choice method for covariance function is GCV if bwxcov = c(0,0)
#' @param ntest1 : number of curves used for CV when choosing bandwidth 
#' @param ngrid1 : number of support points for the covariance surface 
#' @param selection_k : the method of choosing the number of principal components K
#' @param FVE_threshold : Fraction-of-Variance-Explained
#' @param maxk : maximum number of principal components to consider
#' @param regular : do we have regular or sparse functional data
#' @param error : error assumption with measurement error
#' @param ngrid : number of support points in each direction of covariance surface 
#' @param method : method to estimate the PC scores
#' @param shrink : apply shrinkage to estimates of random coefficients (regular data only)
#' @param newdata : new data points to estimate
#' @param kernel : smoothing kernel choice
#' @param numBins : number of bins
#' @param yname : name of the variable analysed
#' @param screePlot : make scree plot
#' @param designPlot : make design plot
#' @param corrPlot : make correlation plot
#' @param rho : truncation threshold for the iterative residual  
#' @param verbose : display diagnostic messages
#' @param xmu : user-defined smoothed mean function 
#' @param xcov : user-defined smoothed covariance function
#' @param method_mu :  method to estimate mu
#' @param out_percent : number in [0,1] indicating the out_percent data in the boundary
#' @param use_binned_data : 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin)
#' @return an option list
#' @examples 
#' 1 + 3

SetOptions = function(y, t, p){ 

  bwmu =p[['bwmu']];                bwmu_gcv =p[['bwmu_gcv']]; 
  bwxcov =p[['bwxcov']];            bwxcov_gcv =p[['bwxcov_gcv']];
  ntest1 =p[['ntest1']];            ngrid1 =p[['ngrid1']]; 
  selection_k =p[['selection_k']];  FVE_threshold =p[['FVE_threshold']];
  maxk =p[['maxk']];                
  regular =p[['regular']];          error =p[['error']];
  ngrid =p[['ngrid']];              method =p[['method']];
  shrink =p[['shrink']];            newdata =p[['newdata']] ;
  kernel =p[['kernel']];            numBins =p[['numBins']];
  yname =p[['yname']];              screePlot =p[['screePlot']];
  designPlot =p[['designPlot']];    rho =p[['rho']];
  verbose =p[['verbose']];          corrPlot =p[['corrPlot']];
  xmu =p[['xmu']];                  method_mu =p[['method_mu']];
  out_percent =p[['out_percent']];  xcov =p[['xcov']];
                                    use_binned_data =p[['use_binned_data']];

  if(is.null(bwmu)){ # bandwidth choice for mean function is using CV or GCV
    bwmu = 0;   
  }
  if(is.null(bwmu_gcv)){ # bandwidth choice for mean function is GCV if bwmu = 0
    bwmu_gcv = 'GMeanAndGCV';  
  }
  if(is.null(bwxcov)){ # bandwidth choice for covariance function is CV or GCV
    bwxcov = 0; 
  }
  if(is.null(bwxcov_gcv)){  # bandwidth choice for covariance function is GCV if bwxcov = c(0,0)
    bwxcov_gcv = 'GMeanAndGCV';
  }
  if(is.null(ntest1)){ # number of curves used for CV when choosing bandwidth 
    ntest1 = min(30, length(y)-1);
  }
  if(is.null(ngrid1)){ # number of support points for the covariance surface 
    ngrid1 = 30;
  }
  if(is.null(selection_k)){ # the method of choosing the number of principal components K
    selection_k = "AIC";
  }
  if(selection_k == "FVE" && is.null(FVE_threshold)){  # the Fraction-of-Variance-Explained
    FVE_threshold = 0.95;
  }
  if(is.null(maxk)){ # maximum number of principal components to consider
    maxk = min(20, length(y)-1);   
  }
  if(is.null(regular)){ #do we have regular or sparse functional data
    regular = IsRegular(t);    
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
    if(regular == "Dense"){
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
  # if(error == FALSE && (selection_k == "AIC" || selection_k == "BIC")){ # Check suitability of information criterion
  #  cat('When assume no measurement error, cannot use "AIC" or "BIC". Reset to "BIC" now!\n')
  #  selection_k = "BIC" 
  #}
  kernNames = c("rect", "gauss", "epan", "gausvar", "quar");
  if(!(kernel %in% kernNames)){ # Check suitability of kernel
    cat(paste('kernel', kernel, 'is unrecognizable! Reset to Epanechnikov kernel now!\n')); 
    kernel = "epan"; 
  }  
  if(!is.null(numBins)){ 
    if(numBins < 10){   # Check suitability of number of bins
      cat("Number of bins must be at least +10!!\n");
      numBins = NULL;
    }
  }
  if(is.null(use_binned_data)){ 
    use_binned_data = 'AUTO';
  }
    
  return( list(bwmu = bwmu, bwmu_gcv = bwmu_gcv, bwxcov = bwxcov, bwxcov_gcv = bwxcov_gcv,
          ntest1 = ntest1, ngrid1 = ngrid1, selection_k = selection_k, FVE_threshold = FVE_threshold,
          maxk = maxk, regular = regular, error = error, ngrid = ngrid, 
          method = method, shrink = shrink, newdata = newdata, kernel = kernel, corrPlot = corrPlot,	
          numBins = numBins, yname = yname, screePlot = screePlot, designPlot = designPlot, rho = rho, 
          verbose = verbose, xmu = xmu, xcov = xcov, method_mu= method_mu, out_percent = out_percent, use_binned_data = use_binned_data) )
}
