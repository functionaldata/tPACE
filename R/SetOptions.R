#' Set the PCA option list
#'
#' @param y A list of \emph{n} vectors containing the observed values for each individual.
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' See '?FPCA for more details. Casual users are not advised to tamper with this function.
#' @export


SetOptions = function(y, t, optns){

  methodMuCovEst = optns[['methodMuCovEst']]
  userBwMu =optns[['userBwMu']];                
  methodBwMu =optns[['methodBwMu']]; 
  userBwCov =optns[['userBwCov']];            
  methodBwCov =optns[['methodBwCov']];
  kFoldMuCov = optns[['kFoldMuCov']]
  methodSelectK =optns[['methodSelectK']];  
  FVEthreshold =optns[['FVEthreshold']];
  fitEigenValues <- optns[['fitEigenValues']];
  maxK =optns[['maxK']];                
  dataType =optns[['dataType']];          
  error =optns[['error']];
  nRegGrid =optns[['nRegGrid']];              
  methodXi =optns[['methodXi']];
  shrink =optns[['shrink']]
  kernel =optns[['kernel']];            
  numBins =optns[['numBins']];
  yname =optns[['yname']];
  rho =optns[['rho']];
  usergrid =optns[['usergrid']];
  userRho = optns[['userRho']];
  diagnosticsPlot =optns[['diagnosticsPlot']];
  plot =optns[['plot']]
  if (!is.null(diagnosticsPlot)) {
    warning("The option 'diagnosticsPlot' is deprecated. Use 'plot' instead")
    plot = diagnosticsPlot
  } 
  verbose =optns[['verbose']];   
  userMu =optns[['userMu']];                  
  #methodMu =optns[['methodMu']];
  outPercent =optns[['outPercent']];  
  userCov =optns[['userCov']];
  userSigma2 = optns[['userSigma2']]
  rotationCut =optns[['rotationCut']];    
  useBinnedData =optns[['useBinnedData']];
  useBinnedCov = optns[['useBinnedCov']]
  lean = optns[['lean']]
  useBW1SE =optns[['useBW1SE']];   

  if(is.null(methodBwMu)){ # bandwidth choice for mean function is GCV if userBwMu = 0
    #methodBwMu = 'GMeanAndGCV';  
    methodBwMu = 'Default'
  }
  if(is.null(userBwMu) && methodBwMu == 'Default'){ # bandwidth choice for mean function is using CV or GCV
    userBwMu = 0.05 * diff(range(unlist(t)));   
  } 
  if(is.null(userBwMu) && methodBwMu != 'Default'){
    userBwMu = 0.0;
  }
  if(is.null(methodBwCov)){  # bandwidth choice for covariance function is GCV if userBwCov = c(0,0)
    #methodBwCov = 'GMeanAndGCV';
    methodBwCov = 'Default';
  }
  if(is.null(userBwCov) && methodBwCov == 'Default'){ # bandwidth choice for covariance function is CV or GCV
    userBwCov = 0.10 * diff(range(unlist(t))); 
  }
  if(is.null(userBwCov) && methodBwCov != 'Default'){
    userBwCov = 0.0;
  }
  #if(is.null(ngrid1)){ # number of support points for the covariance surface 
  #  ngrid1 = 30;
  #}
  if (is.null(kFoldMuCov)) { # CV fold for covariance smoothing
    kFoldMuCov <- 10L
  } else {
    kFoldMuCov <- as.integer(kFoldMuCov)
  }
  if(is.null(methodSelectK)){ # the method of choosing the number of principal components K
    #  TODO : Possibly have user-defined selection methods for the # of FPCs and we keep
    # an internal FVE-based method for us
    methodSelectK = "FVE";
  }
  if(is.null(FVEthreshold)){  # Default Value for the Fraction-of-Variance-Explained
     FVEthreshold = 0.9999;
  }
  if(is.null(dataType)){ #do we have dataType or sparse functional data
    dataType = IsRegular(t);    
  }
  if (is.null(fitEigenValues)) {
    fitEigenValues <- FALSE
  }
  if(is.null(methodMuCovEst)){
    if (dataType == 'Sparse'){
      methodMuCovEst = 'smooth'; #In the odd case that somehow we use this...
    } else {
      methodMuCovEst = 'cross-sectional';
    }
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
  if(is.null(maxK)){ # maximum number of principal components to consider
    maxK = min( nRegGrid-2, length(y)-2 );   
    if(methodMuCovEst == 'smooth'){
      maxK = min( maxK, 20) 
    }
    if(maxK < 1){
      message("Automatically defined maxK cannot be less than 1. Reset to maxK = 1 now!\n")
      maxK = 1
    }
    if( length(y) <= 3 ){
      message("The sample size is less or equal to 3 curves. Be cautious!\n")
    }
  }
  methodNames = c("IN", "CE");
  if(!is.null(methodXi) && !(methodXi %in% methodNames)){
    message(paste('methodXi', methodXi, 'is unrecognizable! Reset to automatic selection now!\n')); 
    methodXi = NULL; 
  }   
  if(is.null(methodXi)){ # method to estimate the PC scores
    if(dataType == 'Dense'){
      methodXi = "IN";
    } else if(dataType == 'Sparse'){
      methodXi = "CE";
    } else if(dataType == 'DenseWithMV'){
      methodXi = "CE"; # We will see how IN can work here
    } else { # for dataType = p>>n
      methodXi = "IN";
    }
  }
   if(is.null(shrink)){ 
     # apply shrinkage to estimates of random coefficients (dataType data
     # only)
     shrink = FALSE;
   }
   if(shrink == TRUE && (error != TRUE || methodXi != "IN")){ 
     # Check for valid shrinkage choice
     message('shrinkage method only has effects when methodXi = "IN" and error = TRUE! Reset to shrink = FALSE now!\n');
     shrink = FALSE      
   }
  if(is.null(kernel)){ # smoothing kernel choice
    if(dataType == "Dense"){
      kernel = "epan";   # kernel: Epanechnikov
    }else{
      kernel = "gauss";  # kernel: Gaussian
    }
  }
  kernNames = c("rect", "gauss", "epan", "gausvar", "quar");
  if(!(kernel %in% kernNames)){ # Check suitability of kernel
    message(paste('kernel', kernel, 'is unrecognizable! Reset to automatic selection now!\n')); 
    kernel = NULL; 
  }  
  if(is.null(kernel)){ # smoothing kernel choice
    if(dataType %in% c( "Dense", "DenseWithMV")){
      kernel = "epan";   # kernel: Epanechnikov
    }else{
      kernel = "gauss";  # kernel: Gaussian
    }
  }
  if(is.null(yname)){ # name of the variable analysed
    yname = as.character(substitute(y))      
  }
  if(maxK > (nRegGrid-2)){ # check if a reasonable number of eigenfunctions is requested
    message(paste("maxK can only be less than or equal to", nRegGrid-2,"! Reset to be", nRegGrid-2, "now!\n"));
    maxK = nRegGrid -2;
  }
  if(is.numeric(methodSelectK)){
    FVEthreshold <- 1 # disable FVE selection.
    if(methodSelectK > (nRegGrid-2)){ # check if a reasonable number of eigenfunctions is requested
      message(paste("maxK can only be less than or equal to", nRegGrid-2,"! Reset to be", nRegGrid-2, "now!\n"));
      maxK = nRegGrid -2;
    }else if(methodSelectK <= 0){ # check if a positive number of eigenfunctions is requested
      message("methodSelectK must be a positive integer! Reset to BIC now!\n");
      methodSelectK = "BIC"
      FVEthreshold = 0.95;
    }
  }
  if(is.null(plot)){ # make corrplot
    plot = FALSE;
  }
  if(is.null(rho)){ # truncation threshold for the iterative residual that is used
    # no regularization if sigma2 is specified or assume no measurement error.
    if (!is.null(userSigma2) || error == FALSE) { 
      rho <- 'no'
    } else {
      rho <- 'cv'
    }
  }
  if(is.null(userRho)){
    userRho = NULL
  }
  if(is.null(verbose)){ # display diagnostic messages
    verbose = FALSE;
  }  
  if(is.null(userMu)){ # user-defined mean functions valued at distinct input time points
    userMu <- NULL
  }
  if(is.null(userCov)){
    userCov <- NULL
  }
  if(is.null(outPercent)){ 
    outPercent <- c(0,1)
  }  
  if(is.null(rotationCut)){ 
    rotationCut <- c(0.25,.75)
  } 
  # if(error == FALSE && (methodSelectK == "AIC" || methodSelectK == "BIC")){ # Check suitability of information criterion
  #  message(paste0('When assume no measurement error, cannot use "AIC" or "BIC". Reset to "BIC" now!\n'))
  #  methodSelectK = "BIC" 
  #}
  if(!is.null(numBins)){ 
    if(numBins < 10){   # Check suitability of number of bins
      message("Number of bins must be at least +10!!\n");
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
  if(is.null(usergrid)){ 
    usergrid = TRUE;
  }
  if(is.null(lean)){ 
    lean = FALSE;
  }
  if(is.null(useBW1SE)){ 
    useBW1SE = FALSE;
  }
  # if (!all.equal(outPercent, c(0, 1)) && methodMuCovEst == 'cross-sectional') {
    # stop('outPercent not supported for cross-sectional covariance estimate')
  # }
    
  retOptns <- list(userBwMu = userBwMu, methodBwMu = methodBwMu, userBwCov = userBwCov, methodBwCov = methodBwCov,
          kFoldMuCov = kFoldMuCov, methodSelectK = methodSelectK, FVEthreshold = FVEthreshold,
          fitEigenValues = fitEigenValues, maxK = maxK, dataType = dataType, error = error, shrink = shrink,
          nRegGrid = nRegGrid, rotationCut = rotationCut, methodXi = methodXi, kernel = kernel, 
          lean = lean, diagnosticsPlot = diagnosticsPlot, plot=plot, numBins = numBins, useBinnedCov = useBinnedCov, 
          usergrid = usergrid, yname = yname,  rho = rho, verbose = verbose, userMu = userMu, userCov = userCov, methodMuCovEst = methodMuCovEst,
          userRho = userRho, userSigma2 = userSigma2, outPercent = outPercent, useBinnedData = useBinnedData, useBW1SE = useBW1SE)

  invalidNames <- !names(optns) %in% names(retOptns)
  if (any(invalidNames)) {
    stop(sprintf('Invalid option names: %s',
                 paste0(names(optns)[invalidNames], collapse=', ')))
  }
  return( retOptns )
}
