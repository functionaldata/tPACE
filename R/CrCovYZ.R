# Calculate the raw and the smoothed cross-covariance between functional
# and scalar predictors using bandwidth bw or estimate that bw using GCV
# Ly      : List of N vectors with amplitude information
# Lt      : List of N vectors with timing information
# Ymu     : vector Q-1 
# bw      : scalar
# Z       : vector N-1
# Zmu     : scalar
# returns : list with: 1. smoothed cross-covariance (vector M-1), 2. raw
# cross-covariance (vector N-1) and the bandwidth used for smoothing (scalar)
CrCovYZ <- function(bw = NULL, Z, Zmu = NULL, Ly, Lt = NULL, Ymu = NULL, support = NULL){
  
  # If only Ly and Z are available assume DENSE data
  if( is.matrix(Ly) && is.null(Lt) && is.null(Ymu) ){
    rawCC <- GetRawCrCovFuncScal(Ly = Ly, Z = Z)
    return ( list(smoothedCC = NULL, rawCC = rawCC, bw = bw, score = NULL) )  
  }
  # Otherwise assume you have SPARSE data
  if( is.null(Zmu) ){
    Zmu = mean(Z);
  }  
  # Get the Raw Cross-covariance 
  ulLt = unlist(Lt)
  if (is.null(support) ){
    obsGrid = sort(unique(ulLt))
  } else {
    obsGrid = support
  }
  rawCC = GetRawCrCovFuncScal(Ly = Ly, Lt = Lt, Ymu = Ymu, Z = Z, Zmu = Zmu )

  # If the bandwidth is known already smooth the raw CrCov
  if( is.numeric(bw) ){
    smoothedCC <- smoothRCC(rawCC, bw, obsGrid )
    score = GCVgauss1D( smoothedY = smoothedCC, smoothedX = obsGrid, 
                        rawX = rawCC$tpairn, rawY = rawCC$rawCCov, bw = bw)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = bw, score = score) )
  # If the bandwidth is unknown use GCV to take find it
  } else {
    # Construct candidate bw's
    h0 = 2.0 * minb( sort(ulLt), 2+1); # 2x the bandwidth needed for at least 3 points in a window
    r = diff(range(ulLt))    
    q = (r/(4*h0))^(1/9);   
    bwCandidates = sort(q^(0:19)*h0);
    # Find their associated GCV scores
    gcvScores = rep(Inf, length(bwCandidates))
    for (i in 1:length(bwCandidates)){
      smoothedCC <- try(silent=TRUE, smoothRCC(rawCC, bw = bwCandidates[i], xout = obsGrid ))
      if( is.numeric(smoothedCC) ){
        gcvScores[i] = GCVgauss1D( smoothedY = smoothedCC, smoothedX = obsGrid, 
                                   rawX = rawCC$tpairn, rawY = rawCC$rawCCov, bw = bwCandidates[i])
      }
    }
#browser()
    # Pick the one with the smallest score
    bInd = which(gcvScores == min(gcvScores, na.rm=TRUE));
    bOpt = max(bwCandidates[bInd]);
    smoothedCC <- smoothRCC( rawCC, bw = bOpt, obsGrid )
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = bOpt, score = min(gcvScores, na.rm=TRUE)) )
  }  
}

# Calculate the smooth Covariances between functional and scalar predictors
# rCC     : raw cross covariance list object returned by GetRawCrCovFuncScal
# bw      : scalar
# xout    : vector M-1
# returns : vector M-1
smoothRCC <- function(rCC,bw,xout){
  x = matrix( unlist(rCC),  ncol=2)
  x= x[order(x[,1]),]
  return( tPACE::Rlwls1d(bw=bw, win=rep(1,nrow(x)), yin=x[,2], xin=x[,1], 'gauss', xout=xout) ) 
}

# Calculate GCV cost off smoothed sample assuming a Gaussian kernel
# smoothedY : vector M-1 
# smoothedX : vector M-1
# rawX      : vector N-1
# rawY      : vector N-1
# bw        : scalar
# returns   : scalar
GCVgauss1D <- function( smoothedY, smoothedX, rawX, rawY, bw){
  cvsum = sum( (rawY - approx(x=smoothedX, y=smoothedY, xout=rawX)$y)^2 );
  k0 = 0.398942;  # hard-coded constant for Gaussian kernel
  N  = length(rawX)
  r  = diff(range(rawX)) 
  return( cvsum / (1-(r*k0)/(N*bw))^2 )
}
