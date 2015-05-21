gcvlwls1d0 <- function(yy,tt, kernel, npoly, nder, regular, verbose=TRUE) {
# This function computes the optimal bandwidth choice for the mean
# function use GCV method by pooling the longitudinal data together. 
# verbose is unused for now
# this is incompatible with PACE because the GCV is calculated in a different way
  
  t = unlist(tt);
  y = unlist(yy);
  
  r = diff(range(t))
  N = length(t);
  
  # Specify the starting bandwidth candidates
  if ( regular == "Sparse") {
    dstar = minb(t, npoly+2);
    if ( dstar > r*0.25){
      dstar = dstar * 0.75;
      warning( c( "The min bandwidth choice is too big, reduce to ", dstar, "!\n"))
    }
    h0 = 2.5 * dstar;
  }else if(regular == "DenseWithMV"){
    h0 = 2.0 * minb(t, npoly+1);
  } else {
    h0 = 1.5 * minb(t,npoly+1);
  }  
  if ( is.nan(h0) ){
    if ( kernel == "gauss" ){
      h0 = 0.2 * r;
    }else{
      error("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
    }
  }  
  h0 <- min(h0,r)
  q = (r/(4*h0))^(1/9); 
  bwCandidates = sort(q^(0:9)*h0) ;
 
  # Get the corresponding GCV scores 
  gcvScores <- unlist(lapply(1:length(bwCandidates), function(i) 
                       gcv(lwls1d(bwCandidates[i], kern=kernel, npoly=npoly, nder=nder, xin = t, yin= y, returnFit=TRUE))[4]))  
  
  # If no bandwith gives a finite gcvScore increase the candidate bandwith and retry on a finer grid
  if(all((is.infinite(gcvScores)))){
    bwCandidates = seq( max(bwCandidates), r, length.out = 2*length(bwCandidates))
    gcvScores <- unlist(lapply(1:length(bwCandidates), function(i)
                        gcv(lwls1d(bwCandidates[i], kern=kernel, npoly=npoly, nder=nder, xin = t, yin= y, returnFit=TRUE))[4]))
  }

  # If the problem persist we clearly have too sparse data
  if(all((is.infinite(gcvScores)))){
    error("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")      
  }

  bInd = which(gcvScores == min(gcvScores));
  bScr = gcvScores[bInd]
  bOpt = bwCandidates[bInd]; 
  
  if( bOpt == r){
    warning("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
  }
  bOptList <- list( 'bOpt' = bOpt, 'bScore' = bScr) 
  return( bOptList)
}
