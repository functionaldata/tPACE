#' Minimum bandwidth based on kNN criterion.
#'
#' Input a list of time points Lt, and the number of unique neighbors k. Obtain  the minimum bandwidth  guaranteeing k unique neighbours.
#' 
#' @param Lt n-by-1 list of vectors  
#' @param k number of unique neighbors for cov and mu (default = 3)
#' @param onlyCov Indicator to return only the minimum bandwidth for the covariance
#' @param onlyMean Indicator to return only the minimum bandwidth for the mean
#' @examples 
#' tinyGrid = list(c(1,7), c(2,3),  6,  c(2,4), c(4,5))
#' BwNN(tinyGrid, k = 2) # c(3,2)
#' @export

BwNN <- function(Lt, k=3, onlyMean = FALSE, onlyCov = FALSE) {
  
  tPairs <- do.call(rbind, lapply(Lt, function(t) {
    expand.grid(t, t)
  }))
  
  if( k <1){
    stop("You cannot have less than 1 neighbours.")
  } 
  
  if( onlyMean && onlyCov){
    stop("BwNN returns nothing!")
  }
  
  distNN2 = NULL
  distNN1 = NULL
  
  if( !onlyMean ){
    uniqTPairs <- unique(tPairs)
    distNN2 <- FindNN(uniqTPairs, k)
  }
  
  if( !is.null(distNN2) && is.infinite(distNN2)){
    stop("You are asking an unreasonable ammount of neighbours for the covariace.")
  }
  
  if( !onlyCov ){
    gridPts <- sort(unique(uniqTPairs[, 1]))
    distNN1 <- max(diff(gridPts, lag=k))
  }
  
  if( !is.null(distNN1) && is.infinite(distNN1)){
    stop("You are asking an unreasonable ammount of neighbours for the mean.")
  }
  
  return(c(cov = distNN2, mu = distNN1))
}

FindNN <- function(mat, k=3) {
  max(apply(mat, 1, function(x) {
    d <- abs(mat - matrix(1, nrow(mat), 1) %*% x)
    distTox <- pmax(d[, 1], d[, 2])
    sort(distTox, partial=k + 1)[k + 1]
  }))
}

