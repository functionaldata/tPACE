#' Sparsify densely observed functional data
#'
#' Given a matrix of densely observed functional data, make a sparsified sample.
#' 
#' @param samp A matrix of densely observed functional data, with each row containing one sample.
#' @param pts A vector of grid points corresponding to the columns of \code{samp}.
#' @param sparsity A vector of integers. The number of observation per sample is chosen to be one of the elements in sparsity with equal chance.
#' @param aggressive Sparsify in an "aggressive" manner making sure that near-by readings are excluded.
#'
#' @return A list of length 2, containing the following fields:
#' \item{tList}{A list of observation time points for each sample.}
#' \item{yList}{A list of values for each sample, corresponding to the time points.}
#' @export
sparsify <- function(samp, pts, sparsity, aggressive = FALSE) {
    if (length(sparsity) == 1)
            sparsity <- c(sparsity, sparsity) # avoid scaler case
    if(!aggressive){
      indEach <- lapply(1:nrow(samp), function(x) 
        sort(sample(ncol(samp), sample(sparsity, 1))))
    } else {
      indEach <- lapply(1:nrow(samp), function(x)
        remotesampling(ncol(samp), sparsity) )
    }
    tList <- lapply(indEach, function(x) pts[x])
    yList <- lapply(1:length(indEach), function(x) {
        ind <- indEach[[x]]
        y <- samp[x, ind]
        return(y)
    })
   
    return(list(tList=tList, yList=yList))
}



remotesampling <- function(N,s){
  onesamp = sort( sample( N, sample( s, 1))) 
  threshold = (1/length(onesamp))^(1.5) * N 
  while( min(diff(onesamp)) < threshold ){
    onesamp = sort( sample( N, length(onesamp)))
  }
  return(onesamp)
}



