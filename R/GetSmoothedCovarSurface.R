GetSmoothedCovarSurface <- function(y, t, out1, mu, p){

  regular = p$regular
  error = p$error
   
  # Get raw covariance   
  rcov <- rcov <-GetRawCov(y, t, out1, mu, regular, error)

  if (any(p$bwxcov ==0)){
    ... 
    }
   
  }

