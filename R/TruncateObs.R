TruncateObs <- function(y, t, obsGrid, buff=.Machine$double.eps * max(abs(obsGrid)) * 3) {

  tmpInd <- mapply(function(yVec, tVec) {
                  ind <- (tVec >= min(obsGrid) - buff & tVec <= max(obsGrid) + buff)
                  return(ind)
                }, y, t)
  y <- mapply(function(yVec, ind) yVec[ind], y, tmpInd)
  t <- mapply(function(tVec, ind) tVec[ind], t, tmpInd)

  return(list(y=y, t=t))
}
