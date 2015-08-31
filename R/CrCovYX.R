# Calculate the raw and the smoothed cross-covariance between functional
# predictors using bandwidths bw's or estimate these bw's using GCV
# Ly1      : List of N vectors with amplitude information
# Lt1      : List of N vectors with timing information
# Ymu1     : vector Q-1 
# bw1      : scalar
# Ly2      : List of N vectors with amplitude information
# Lt2      : List of N vectors with timing information
# Ymu2     : vector Q-1 
# bw2      : scalar
# returns : list with: 1. smoothed cross-covariance (Matrix M-M), 2. raw
# cross-covariance (vector N-1) and 3. the bandwidths used for smoothing (vectors 2-1)
CrCovYX <- function(bw1 = NULL, bw2 = NULL, Ly1, Lt1 = NULL, Ymu1 = NULL, Ly2, Lt2 = NULL, Ymu2 = NULL){
  
  # If only Ly1 and Ly2 are available assume DENSE data
  if( is.matrix(Ly1) && is.null(Lt1) && is.null(Ymu1) && is.matrix(Ly2) && is.null(Lt2) && is.null(Ymu2)){
    rawCC <- GetRawCrCovFuncFunc(Ly1 = Ly1, Ly2 = Ly2)
    return ( list(smoothedCC = NULL, rawCC = rawCC, bw = NULL, score = NULL) )  
  }
  
  # Otherwise assume you have SPARSE data
  if( is.null(Ymu1) ||   is.null(Ymu2)){
    stop("Both functional means must be provided.")   
  }  
  
  # Get the Raw Cross-covariance    
  rawCC = GetRawCrCovFuncFunc(Ly1 = Ly1, Lt1 = Lt1, Ymu1 = Ymu1, Ly2 = Ly2, Lt2 = Lt2, Ymu2 = Ymu2)
  
  # Calculate the observation and the working grids
  ulLt1 = unlist(Lt1);             ulLt2 = unlist(Lt2)
  obsGrid1 = sort(unique(ulLt1));  obsGrid2 = sort(unique(ulLt2))
  
  workGrid1 = seq(obsGrid1[1], max(obsGrid1), length.out = 50)
  workGrid2 = seq(obsGrid2[1], max(obsGrid2), length.out = 50)
  workGrid12 = matrix(c(workGrid1, workGrid2),ncol= 2)
  # If the bandwidth is known already smooth the raw CrCov
  if( is.numeric(bw1) &&  is.numeric(bw2)){
    smoothedCC <- smoothRCC2D(rcov =rawCC, bw1, bw2, workGrid1, workGrid2)    
    score = GCVgauss2D(smoothedCC = smoothedCC, smoothGrid = workGrid12, 
                       rawCC = rawCC$rawCCov, rawGrid = rawCC$tpairn, bw1 = bw1, bw2 = bw2)                      
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw =  c(bw1, bw2), score = score, smoothGrid = workGrid12 ) )
    
  # If the bandwidths are unknown use GCV to take find it
  } else {
    # Construct candidate bw's
    bwCandidates <- getBWidths(ulLt1, ulLt2)
    # Find their associated GCV scores 
    gcvScores = rep(Inf, nrow(bwCandidates)) 
    for (i in 1:length(bwCandidates)){
      smoothedCC <- try(silent=TRUE, smoothRCC2D(rcov=rawCC, bw1 = bwCandidates[i,1], 
                                                 bw2 = bwCandidates[i,2], workGrid1, workGrid2) )
      if( is.numeric(smoothedCC) ){
        gcvScores[i] = GCVgauss2D( smoothedCC = smoothedCC, smoothGrid = workGrid12, rawCC = rawCC$rawCCov, 
                       rawGrid = rawCC$tpairn, bw1 = bwCandidates[i,1], bw2 = bwCandidates[i,2])
      }
    } 
    # Pick the one with the smallest score
    bInd = which(gcvScores == min(gcvScores, na.rm=TRUE));
    bOpt1 = max(bwCandidates[bInd,1]);
    bOpt2 = max(bwCandidates[bInd,2]); 
    smoothedCC <- smoothRCC2D(rcov=rawCC, bw1 =bOpt1, bw2 =bOpt2, workGrid1, workGrid2)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = c(bOpt1, bOpt2), smoothGrid = workGrid12, score = min(gcvScores, na.rm=TRUE)) )
  }  
}

getBWidths <- function(ulLt1, ulLt2){

  bwCandidates <- matrix(rep(0,72),ncol=2)
  h0 = 2.0 * minb( sort(ulLt1), 2+1); # 2x the bandwidth needed for at least 3 points in a window
  r = diff(range(ulLt1))    
  q = (r/(4*h0))^(1/9);  
  bwCandidates[,1] = rep( sort(q^( seq(0,12,length.out=6) )*h0), times= 6);
  h0 = 2.0 * minb( sort(ulLt2), 2+1); # 2x the bandwidth needed for at least 3 points in a window
  r1 = diff(range(ulLt2))    
  q = (r/(4*h0))^(1/9);  
  bwCandidates[,2] =  rep( sort(q^( seq(0,12,length.out=6) )*h0), each= 6); 
  
  return(bwCandidates)
}

# Calculate the smooth Covariances between two functional variables
# rcov    : raw cross covariance list object returned by GetRawCrCovFuncFunc
# bw1     : scalar
# bw2     : scalar
# xout1   : vector M-1
# xout2   : vector L-1
# returns : matrix M-L
smoothRCC2D <- function(rcov,bw1, bw2, xout1, xout2){
  return( tPACE::lwls2dV2( bw = c(bw1, bw2), kern = 'gauss', xin=rcov$tpairn, 
                           yin=rcov$rawCC, xout1=xout1, xout2=xout2, crosscov=TRUE) )  
}

# Calculate GCV cost off smoothed sample assuming a Gaussian kernel
# smoothedY : vector M-1 
# smoothedX : vector M-1
# rawX      : vector N-1
# rawY      : vector N-1
# bw        : scalar
# returns   : scalar
GCVgauss2D <- function( smoothedCC, smoothGrid, rawCC, rawGrid, bw1, bw2){ 

  obsFit <- tPACE::interp2lin(smoothGrid[,1], smoothGrid[,2], smoothedCC, rawGrid[, 1], rawGrid[, 2])    
  # workaround for degenerate case.
  if (any(is.nan(obsFit))  || any(is.infinite(obsFit))  ){
    return(Inf)
  }  
  # residual sum of squares
  cvsum <- sum((rawCC - obsFit) ^ 2)
  N   = length( rawCC )
  r1  = diff( range(smoothGrid[,1] ) )
  r2  = diff( range(smoothGrid[,2] ) )
  k0 = 0.398942;  # hard-coded constant for Gaussian kernel
  return( cvsum / (1 - (1/N) * (r1 * k0 * r2 * k0) /(bw1 * bw2))^2 )
}
