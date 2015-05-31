BinRawCov <- function(rcov) {
  # rcov: Assumes rcov has entries on a regular grid.
  # Returns: 
  
  rcov$tpairn <- round(rcov$tpairn, 13)
  
  # Get the count, mean raw cov, and residual sum of squares at each pair of observed time points.
  tmp <- tapply(rcov$cxxn, list(rcov$tpairn[, 1], rcov$tpairn[, 2]), function(yy) c(mean(yy), length(yy), var(yy) * (length(yy) - 1)), simplify=FALSE)
  # Corresponding time pairs to the non null values
  tPairs <- unname(as.matrix(expand.grid(as.numeric(dimnames(tmp)[[1]]), as.numeric(dimnames(tmp)[[2]]))[!sapply(tmp, is.null), ]))
  tmp <- do.call(rbind, tmp)
  meanVals <- tmp[, 1]
  count <- tmp[, 2]
  RSS <- tmp[, 3] # Residual sum of squares. For implementing GCV.
  RSS[is.na(RSS)] <- 0
  
  diagRSS <- diagCount <- diagMeans <- tDiag <- NULL
  if (!is.null(rcov$diag)) {
    tmp <- do.call(rbind, tapply(rcov$diag[, 2], rcov$diag[, 1], function(yy) c(mean(yy), length(yy), var(yy) * (length(yy) - 1))))
    tDiag <- as.numeric(rownames(tmp))
    diagMeans <- unname(tmp[, 1])
    diagCount <- unname(tmp[, 2])
    diagRSS <- unname(tmp[, 3])
    diagRSS[is.na(diagRSS)] <- 0
  }
  
  res <- list(tPairs=tPairs, meanVals=meanVals, RSS=RSS, tDiag=tDiag, diagMeans=diagMeans, diagRSS=diagRSS, count=count, diagCount=diagCount, error=rcov$error, regular=rcov$regular)
  class(res) <- 'BinnedRawCov'
  
  return(res)
}