### CheckAndCreateCOPoptions # COP = CreateOutliersPlot

CheckAndCreateCOPoptions <- function(optns,fObjClass){
  
  if(is.null(optns$ifactor)){
    ifactor = NULL
  } else {
    ifactor = optns$ifactor
  }
  
  if(is.null(optns$outlierList)){
    outlierList = NULL
  } else {
    outlierList = optns$outlierList
  }
  
  if(is.null(optns$unimodal)){
    unimodal = NULL
  } else {
    unimodal = optns$unimodal
  }
  
  if(is.null(optns$colSpectrum)){
    colSpectrum = NULL
  } else {
    colSpectrum = optns$colSpectrum
  }
  
  if(is.null(optns$groupingType)){
    groupingType = 'standard'
  } else {
    groupingType = optns$groupingType
  }
  
  if(is.null(optns$variant)){
    variant = 'KDE'
  } else {
    variant = optns$variant
  }
  
  if(is.null(optns$nSlices)){
    nSlices = 4
  } else {
    nSlices = optns$nSlices
  }
  
  
  if(is.null(optns$showSlices)){
    showSlices = FALSE
  } else {
    showSlices = optns$showSlices
  }
  
  if(is.null(optns$fIndeces)){
    if(fObjClass == 'FPCA'){
      fIndeces <- c(1,2)
    } else {
      fIndeces <- c(1,1)
    }
  } else {
    if( 2 < length(optns$fIndeces)){
      warning("fIndeces has more than two elements; only the first two will be used.")
    }
    fIndeces <-  (optns$fIndeces[1:2])
  }
  
  if( !any( groupingType == c('standard','slice')) ){
    stop("You request an groupingType method not currenty available.")
  }
  if( !any( variant == c('KDE','bagplot', 'NN')) ){
    stop("You request an outlier detection method not currenty available.")
  }
  if ( variant == 'bagplot' && !is.element('aplpack', installed.packages()[,1]) ){
    stop("Cannot the use the bagplot method; the package 'aplpack' is unavailable.")
  }
  if ( variant == 'KDE' && !is.element('ks', installed.packages()[,1]) ){
    stop("Cannot the use the KDE method; the package 'ks' is unavailable.")
  } 
  if ( !is.null(unimodal) && !is.logical(unimodal) ){
    stop("The variable 'unimodal' must be logical.")
  } 
  if (is.null(colSpectrum)){
    colFunc = colorRampPalette(c("red",  "yellow", 'blue'));
  } else {
    colFunc = colorRampPalette(colSpectrum)
  }
  if (!is.null(ifactor) && (1 >= ifactor) ){
    warning("It is nonsensical for an inflation factor to be <= 1. 'ifactor' set to 1.1.")
    ifactor = 1.1;
  }
  if ( !(2 <= nSlices) || !(16 >= nSlices) || !(nSlices %% 1 == 0) ){
    warning("nSlices must be between a natural number between 2 and 16. 'nSlices' set to 4.")
    nSlices = 4;
  }
  if(diff(range(fIndeces)) < .Machine$double.eps){
    if( fObjClass == 'FPCA'){
      stop("You request a scatter over the same component; check the fIndeces.")
    }
  } 
  if(is.null(optns$maxVar)){ 
    maxVar = !( fObjClass == 'FPCA') 
  } else {
    maxVar = optns$maxVar
  }
  
  
  perfOptns <- list(nSlices = nSlices, ifactor = ifactor, colFunc = colFunc, fIndeces = fIndeces, maxVar = maxVar,
                    showSlices = showSlices,
                    variant = variant, groupingType = groupingType, unimodal = unimodal,  outlierList = outlierList)
  return(perfOptns)
  
  
}