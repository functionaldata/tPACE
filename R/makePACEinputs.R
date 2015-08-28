#' Format FPCA input
#'
#' Turn vector inputs to the list so they can be used in FPCA 
#' 
#' @param IDs  : n-by-1 vector of subject IDs
#' @param tVec : n-by-1 vector of measurement times
#' @param yVec : n-by-1 vector of measurements from the variable of interest
#' @return L   : list containing 3 lists each of length 'm', 'm' being the number of unique subject IDs
#' @export

makePACEinputs <- function(IDs = NULL, tVec, yVec){

  if( !is.null(IDs) ){ 
    uniqueIDs <- unique(IDs) 
    Lt <- lapply( uniqueIDs, function(x) tVec[ which(IDs == x)])
    Ly <- lapply( uniqueIDs, function(x) yVec[ which(IDs == x)])
    Lid <- as.list(uniqueIDs)
  } else if ( is.matrix(yVec) && is.null(IDs) && is.vector(tVec) ){
    Ly <- lapply( seq_len(nrow(yVec)), function(i) yVec[i,])
    Lt <- rep( list(tVec), dim(yVec)[1] )
    Lid <- as.list( 1:dim(yVec)[1])
  }
  L <- list( Lid = Lid, Ly = Ly, Lt = Lt)
  return(L)

}

