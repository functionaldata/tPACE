GetSmoothedCovarSurface <- function(y, t, out1, mu, p){

# TODO: pass in only the options needed, rather than p itself.
  regular = p$regular
  error = p$error
   
  # Get raw covariance   
  rawcov <- GetRawCov(y, t, out1, mu, regular, error)

  if (any(p$bwxcov ==0)){
    ... 
    }
    

  result <- list( rawcov= rawcov, smocov = smocov, bwCov = bwCov );
  class(result) <- "SCS"  
  }

