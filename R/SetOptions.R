#' Set the PCA option list
#'
#' @param y A list of \emph{n} vectors containing the observed values for each individual.
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' See '?FPCAfor more details. Usually users are not supposed to use this function directly.
#'


SetOptions = function(y, t, optns){ 

  bwmu =optns[['bwmu']];                bwmuMethod =optns[['bwmuMethod']]; 
  bwuserCov =optns[['bwcov']];            bwuserCovGcv =optns[['bwuserCovGcv']];
  kFoldCov = optns[['kFoldCov']]
  # ngrid1 =optns[['ngrid1']]; 
  selectionMethod =optns[['selectionMethod']];  FVEthreshold =optns[['FVEthreshold']];
  fitEigenValues <- optns[['fitEigenValues']]
  maxK =optns[['maxK']];                
  dataType =optns[['dataType']];          error =optns[['error']];
  nRegGrid =optns[['nRegGrid']];              methodXi =optns[['methodXi']];
#  shrink =optns[['shrink']];           # newdata =optns[['newdata']] ;
  kernel =optns[['kernel']];            numBins =optns[['numBins']];
#  numComponents =optns[['numComponents']];
  yname =optns[['yname']];             # screePlot =optns[['screePlot']];
#  designPlot =optns[['designPlot']];    
  rho =optns[['rho']];                diagnosticsPlot =optns[['diagnosticsPlot']];
  verbose =optns[['verbose']];          #corrPlot =optns[['corrPlot']];
  userMu =optns[['userMu']];                  methodMu =optns[['methodMu']];
  outPercent =optns[['outPercent']];  userCov =optns[['userCov']];
  rotationCut =optns[['rotationCut']];    useBinnedData =optns[['useBinnedData']];
  #corrPlotType =optns[['corrPlotType']]; 
  useBinnedCov = optns[['useBinnedCov']]

  if(is.null(bwmu)){ # bandwidth choice for mean function is using CV or GCV
    bwmu = 0;   
  }
  if(is.null(bwmuMethod)){ # bandwidth choice for mean function is GCV if bwmu = 0
    bwmuMethod = 'GMeanAndGCV';  
  }
  if(is.null(bwuserCov)){ # bandwidth choice for covariance function is CV or GCV
    bwuserCov = 0; 
  }
  if(is.null(bwuserCovGcv)){  # bandwidth choice for covariance function is GCV if bwuserCov = c(0,0)
    bwuserCovGcv = 'GMeanAndGCV';
  }
  #if(is.null(ngrid1)){ # number of support points for the covariance surface 
  #  ngrid1 = 30;
  #}
  if (is.null(kFoldCov)) { # CV fold for covariance smoothing
    kFoldCov <- 10L
  } else {
    kFoldCov <- as.integer(kFoldCov)
  }
  if(is.null(selectionMethod)){ # the method of choosing the number of principal components K
    #  TODO : Possibly have user-defined selection methods for the # of FPCs and we keep
    # an internal FVE-based method for us
    selectionMethod = "FVE";
  }
  if(selectionMethod == "FVE" && is.null(FVEthreshold)){  # the Fraction-of-Variance-Explained
     FVEthreshold = 0.9999;
  }
  if(is.null(maxK)){ # maximum number of principal components to consider
    maxK = min(20, length(y)-1);   
  }
#  if(is.null(numComponents)){ # maximum number of principal components to return
#    numComponents = NULL;
#  }
  if(is.null(dataType)){ #do we have dataType or sparse functional data
    dataType = IsRegular(t);    
  }
  if (is.null(fitEigenValues)) {
    fitEigenValues <- FALSE
  }
  if (fitEigenValues && dataType == 'Dense') {
    stop('Fit method only apply to sparse data')
  }
  if(is.null(error)){ # error assumption with measurement error
    error = TRUE;    
  }
  if(is.null(nRegGrid)){ # number of support points in each direction of covariance surface 
    if(dataType == 'Dense' || dataType == 'DenseWithMV'){
      tt = unlist(t)
      nRegGrid = length(unique(signif(tt[!is.na(tt)],6)));
    } else { # for Sparse and p>>n
      nRegGrid = 51;
    }    
  }
  methodNames = c("IN", "CE");
  if(!is.null(methodXi) && !(methodXi %in% methodNames)){ # Check suitability of kernel
    cat(paste('methodXi', methodXi, 'is unrecognizable! Reset to automatic selection now!\n')); 
    methodXi = NULL; 
  }   
  if(is.null(methodXi)){ # method to estimate the PC scores
    if(dataType == 'Dense'){
      shrink = FALSE;
      methodXi = "IN";
    } else if(dataType == 'Sparse'){
      shrink = FALSE;
      methodXi = "CE";
    } else if(dataType == 'DenseWithMV'){
      shrink = FALSE;
      methodXi = "IN";
    } else { # for dataType = p>>n
      shrink = FALSE;
      methodXi = "IN";
    }
  }
#  if(is.null(shrink)){ # apply shrinkage to estimates of random coefficients (dataType data only)
#    shrink = FALSE;
#  }
#  if(shrink == TRUE && (error != TRUE || methodXi != "IN")){ # Check for valid shrinkage choice
#    cat('shrinkage method only has effects when methodXi = "IN" and error = TRUE! Reset to shrink = FALSE now!\n');
#    shrink = FALSE      
#  }
  if(is.null(kernel)){ # smoothing kernel choice
    if(dataType == "Dense"){
      kernel = "epan";   # kernel: Epanechnikov
    }else{
      kernel = "gauss";  # kernel: Gaussian
    }
  }
  if(is.null(yname)){ # name of the variable analysed
    yname = as.character(substitute(y))      
  }
  if(maxK > (nRegGrid-2)){ # check if a reasonable number of eigenfunctions is requested
    cat(paste("maxK can only be less than or equal to", nRegGrid-2,"! Reset to be", nRegGrid-2, "now!\n"));
    maxK = nRegGrid -2;
  }
  if(is.numeric(selectionMethod)){
    if(selectionMethod > (nRegGrid-2)){ # check if a reasonable number of eigenfunctions is requested
      cat(paste("maxK can only be less than or equal to", nRegGrid-2,"! Reset to be", nRegGrid-2, "now!\n"));
      maxK = nRegGrid -2;
    }else if(selectionMethod <= 0){ # check if a positive number of eigenfunctions is requested
      cat("selectionMethod must be a positive integer! Reset to BIC now!\n");
      selectionMethod = "BIC"
      FVEthreshold = 0.95;
    }
  }
  #if(is.null(screePlot)){ # make screeplot
  #  screePlot = FALSE;
  #}
  #if(is.null(designPlot)){ # make designplot
  #  designPlot = FALSE;
  #}
  #if(is.null(corrPlot)){ # make corrplot
  #  corrPlot = FALSE;
  #}
  if(is.null(diagnosticsPlot)){ # make corrplot
    diagnosticsPlot = FALSE;
  }
  if(is.null(rho)){ # truncation threshold for the iterative residual that is used
    rho = "cv";
  }
  if(is.null(verbose)){ # display diagnostic messages
    verbose = FALSE;
  }  
  #if(is.null(newdata)){ # new data points to estimate
  #  newdata <- NULL
  #}
  if(is.null(userMu)){ # user-defined mean functions valued at distinct input time points
    userMu <- NULL
  }
  if(is.null(userCov)){
    userCov <- NULL
  }
  if(is.null(methodMu)){ # method to estimate mu
    methodMu <- 'PACE'
  }
  if(is.null(outPercent)){ 
    outPercent <- c(0,1)
  }  
  if(is.null(rotationCut)){ 
    rotationCut <- c(0.25,.75)
  } 
  # if(error == FALSE && (selectionMethod == "AIC" || selectionMethod == "BIC")){ # Check suitability of information criterion
  #  cat('When assume no measurement error, cannot use "AIC" or "BIC". Reset to "BIC" now!\n')
  #  selectionMethod = "BIC" 
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
  if(is.null(useBinnedData)){ 
    useBinnedData = 'AUTO';
  }
  if (is.null(useBinnedCov)) {
    useBinnedCov <- TRUE
    if (  ( 128 > length(y) ) && ( 3 > mean ( unlist( lapply( y, length) ) ) )){
      useBinnedCov <- FALSE
    } 
  }
 
  #if(is.null(corrPlotType)){ 
  #  corrPlotType = 'Fitted';
  #}
    
  return( list(bwmu = bwmu, bwmuMethod = bwmuMethod, bwuserCov = bwuserCov, bwuserCovGcv = bwuserCovGcv,
          kFoldCov = kFoldCov, 
          selectionMethod = selectionMethod, FVEthreshold = FVEthreshold,
          fitEigenValues = fitEigenValues,
          maxK = maxK, dataType = dataType, error = error, nRegGrid = nRegGrid, rotationCut = rotationCut,
          methodXi = methodXi, # shrink = shrink, # newdata = newdata, 
          kernel = kernel, 
          # corrPlot = corrPlot, corrPlotType = corrPlotType,  screePlot = screePlot, designPlot = designPlot,
          # numComponents = numComponents,
          diagnosticsPlot = diagnosticsPlot,
          numBins = numBins, useBinnedCov = useBinnedCov, yname = yname,  rho = rho, 
          verbose = verbose, userMu = userMu, userCov = userCov, methodMu= methodMu, outPercent = outPercent, useBinnedData = useBinnedData) )
}
