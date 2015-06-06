GetSmoothedCovarSurface <- function(y, t, mu, obsGrid, regGrid, optns, useBins=FALSE) {

  # TODO: pass in only the options needed, rather than p itself.
  dataType <- optns$dataType
  error <- optns$error
  kern <- optns$kernel
  bwuserCov <- optns$bwuserCov
  bwuserCovGcv <- optns$bwuserCovGcv
  verbose <- optns$verbose
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  maxGrid <- rangeGrid[2]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] -
                        buff | 
                        regGrid < minGrid + diff(rangeGrid) * outPercent[2] +
                        buff]

  # Get raw covariance   
  rcov <- GetRawCov(y, t, obsGrid, mu, dataType, error)

  if (useBins && bwuserCovGcv == 'CV')
    stop('If bwuserCovGcv == \'CV\' then we must use the unbinned rcov.')

  if (useBins)
    rcov <- BinRawCov(rcov)
  
  if (bwuserCov == 0) { # bandwidth selection
    if (bwuserCovGcv %in% c('GCV', 'GMeanAndGCV')) { # GCV
      gcvObj <- gcvlwls2d(obsGrid, kern=kern, rcov=rcov, verbose=verbose)
      bwCov <- gcvObj$h
      if (bwuserCovGcv == 'GMeanAndGCV') {
        bwCov <- sqrt(bwCov * gcvObj$minBW)
      }  
    } else if (bwuserCovGcv == 'CV') { # CV 10 fold
      gcvObj <- gcvlwls2d(obsGrid, kern=kern, rcov=rcov,
                          verbose=optns$verbose, CV='10fold')
      bwCov <- gcvObj$h
    }
  } else if (bwuserCov != 0) {
    bwCov <- bwuserCov
  }

  if (!useBins)
    smoothCov <- lwls2d(bwCov, kern, xin=rcov$tpairn, yin=rcov$cxxn,
                        xobsGrid=cutRegGrid, xout2=cutRegGrid)
  else 
    smoothCov <- lwls2d(bwCov, kern, xin=rcov$tPairs, yin=rcov$meanVals,
                        win=rcov$count, xobsGrid=cutRegGrid, xout2=cutRegGrid)
  
  # TODO: add cut argument
  if (error)
    sigma2 <- pc_covE(obsGrid, regGrid, bwCov, kernel=kern, rcov=rcov)$sigma2
  else 
    sigma2 <- NULL

  res <- list( rawCov= rcov, smoothCov = (smoothCov + t(smoothCov)) / 2, bwCov = bwCov, sigma2 = sigma2);
  class(res) <- "SmoothCov"  
  
  return(res)
}

