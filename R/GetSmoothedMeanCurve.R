GetSmoothedMeanCurve <- function (y, t, obsGrid, regGrid, optns){
  
  # Note : If binned data we should use weighted mean response for each time-point.
  # This is not currently implemented. \hat{y}_i = \sum_i w_i y_i where w_i are the
  # same points for common t_is. so we have: \hat{y}_i = n_t w_i \bar{y}

  userMu = optns$userMu;
  bwmuGcv = optns$bwmuGcv;
  npoly = 1
  nder = 0 
  bwmu = optns$bwmu; 
  kernel = optns$kernel
 
  # If the user provided a mean function use it
  if (!(isempty(userMu)) && (length(userMu) == length(obsGrid))){
    mu = userMu;
    muDense = interp1(obsGrid,mu, regGrid, 'spline');
    bw_mu = NULL;
  # otherwise if the user provided a mean bandwidth use it to estimate the mean function (below)
  } else {
    if (bwmu > 0){
      bw_mu = bwmu;
    #otherwise estimate the mean bandwith via the method selected to estimnate the mean function (below)
    } else {
      if( any(bwmuGcv == c('GCV','GMeanAndGCV') )){
        # get the bandwidth using GCV
        bw_mu =  unlist(gcvlwls1d1(yy = y, tt = t, kernel = kernel, npoly = npoly, nder = nder, dataType = optns$dataType) )[1]    
        if ( isempty(bw_mu)){ 
          stop('The data is too sparse to estimate a mean function. Get more data!\n')
         }
         bw_mu = adjustBW1(kernel=kernel,bopt=bw_mu,npoly=npoly,dataType=optns$dataType,nder=nder)
         # get the geometric mean between the minimum bandwidth and GCV bandwidth to estimnate the mean function (below)         
         if ( bwmuGcv == 'GMeanAndGCV') {
           minbw = minb( unlist(t),2)
           bw_mu = sqrt(minbw*bw_mu);
        } 
      } else {
        # get the bandwidth using CV to estimnate the mean function (below)
        bw_mu = cvlwls(x=y, tt=t, kernel= kernel, npoly=npoly, nder=nder, dataType= optns$dataType ); 
      }
    }
    # Get the mean function using the bandwith estimated above:
    mu = lwls1d(bw_mu, kern = kernel, npoly = npoly, nder = nder, xin = unlist(t), yin= unlist(y), xout = obsGrid)
    muDense = lwls1d(bw_mu, kern = kernel, npoly = npoly, nder = nder, xin = unlist(t), yin= unlist(y), xout = regGrid)
  }  
  
  result <- list( 'mu' = mu, 'muDense'= muDense, 'bw_mu' = bw_mu);
  class(result) <- "SMC"
  return(result)
}


