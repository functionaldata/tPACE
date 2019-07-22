# TODO: Roxygen documentation

TruncateObs <- function(Ly, Lt, obsGrid, buff=.Machine$double.eps * max(abs(obsGrid)) * 3) { 
  tmpInd <- mapply(function(yVec, tVec) {
                  ind <- (tVec >= min(obsGrid) - buff & tVec <= max(obsGrid) + buff)
                  return(ind)
                }, Ly, Lt, SIMPLIFY=FALSE)
  Ly <- mapply(function(yVec, ind) yVec[ind], Ly, tmpInd, SIMPLIFY = FALSE)
  Lt <- mapply(function(tVec, ind) tVec[ind], Lt, tmpInd, SIMPLIFY = FALSE)

  return(list(Ly=Ly, Lt=Lt))
}
