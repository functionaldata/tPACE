GetSmoothedMeanCurve <- function (y, t, out1, out21, p){
  
  # Note : If binned data we should use weighted mean response for each time-point.
  # This is not currently implemented. \hat{y}_i = \sum_i w_i y_i where w_i are the
  # same points for common t_is. so we have: \hat{y}_i = n_t w_i \bar{y}

  xmu = p$xmu;
  bwmu_gcv = p$bwmu_gcv;
  npoly = 1
  nder = 0 
  bwmu = p$bwmu; 
  kernel = p$kernel
 
  # If the user provided a mean function use it
  if (!(isempty(xmu)) && (length(xmu) == length(out1))){
    mu = xmu;
    muDense = interp1(out1,mu, out21, 'spline');
    bw_mu = NULL;
  # otherwise if the user provided a mean bandwidth use it to estimate the mean function (below)
  } else {
    if (bwmu > 0){
      bw_mu = bwmu;
    #otherwise estimate the mean bandwith via the method selected to estimnate the mean function (below)
    } else {
      if( any(bwmu_gcv == c('GCV','GMeanAndGCV') )){
        # get the bandwidth using GCV
        bw_mu =  gcvlwls1d1(yy = y, tt = t, kernel = kernel, npoly = npoly, nder = nder, regular = p$regular )    
        if ( isempty(bw_mu$bOpt)){ 
          stop('The data is too sparse to estimate a mean function. Get more data!\n')
         }
         # get the geometric mean between the minimum bandwidth and GCV bandwidth to estimnate the mean function (below)
         if ( bwmu_gcv == 'GMeanAndGCV') {
           minbw = minb( unlist(t),2)
           bw_mu = sqrt(minbw*bw_mu$bOpt);
        } 
      } else {
        # get the bandwidth using CV to estimnate the mean function (below)
        bw_mu = cvlwls(x=y, tt=t, kernel= kernel, npoly=npoly, nder=nder, regular= p$regular ); 
      }
    }
    # Get the mean function using the bandwith estimated above:
    mu = lwls1d(bw_mu, kern = kernel, npoly = npoly, nder = nder, xin = unlist(t), yin= unlist(y), xout = out1)
    muDense = lwls1d(bw_mu, kern = kernel, npoly = npoly, nder = nder, xin = unlist(t), yin= unlist(y), xout = out21)
  }  
  
  result <- list( 'mu' = mu, 'muDense'= muDense, 'bw_mu' = bw_mu);
  class(result) <- "SMC"
  return(result)
}

