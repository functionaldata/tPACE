GetSmoothedCovarSurface <- function(y, t, mu, obsGrid, regGrid, p, useBins=FALSE) {

  # TODO: pass in only the options needed, rather than p itself.
  regular <- p$regular
  error <- p$error
  kern <- p$kernel
  bwxcov <- p$bwxcov
  bwxcov_gcv <- p$bwxcov_gcv

  # Get raw covariance   
  rcov <- GetRawCov(y, t, obsGrid, mu, regular, error)

  # TODO: bin rcov
  # If bwxcov_gcv == 'CV' then we must use the unbinned rcov.
  if (useBins && bwxcov_gcv != 'CV')
    rcov <- BinRawCov(rcov)
  
  if (bwxcov == 0) { # bandwidth selection
    if (bwxcov_gcv %in% c('GCV', 'GMeanAndGCV')) { # GCV
      gcvObj <- gcvlwls2d(obsGrid, kern=kern, rcov=rcov, verbose=p$verbose)
      bwCov <- gcvObj$h
      if (bwxcov_gcv == 'GMeanAndGCV') {
        bwCov <- sqrt(bwCov * gcvObj$minBW)
      }  
    } else if (bwxcov_gcv == 'CV') { # CV 10 fold
      gcvObj <- gcvlwls2d(obsGrid, kern=kern, rcov=rcov, verbose=p$verbose, CV='10fold')
      bwCov <- gcvObj$h
    }
  } else if (bwxcov != 0) {
    bwCov <- bwxcov
  }

  smoothCov <- lwls2d(bwCov, kern, xin=rcov$tpairn, yin=rcov$cxxn, xout1=regGrid, xout2=regGrid)
  
  # TODO: add cut argument
  if (error)
    sigma2 <- pc_covE(obsGrid, regGrid, bwCov, kernel=kern, rcov=rcov)$sigma2
  else 
    sigma2 <- NULL

  res <- list( rawCov= rcov, smoothCov = (smoothCov + t(SmoothCov)) / 2, bwCov = bwCov, sigma2 = sigma2);
  class(res) <- "SmoothCov"  
  
  return(res)
}

