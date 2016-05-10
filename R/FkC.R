#' Functional clustering and identifying substructures of longitudinal data
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param k A scalar defining the number of clusters to define; default 3.
#' @param maxIter A scalar defining the maximum number of iterations allowed; default 20, common for both the simple kmeans initially and the subsequent k-centres
#' @param optns A list of options control parameters specified by \code{list(name=value)}; by default: 'maxK' is set to 3. See `Details in ?FPCA'.
#'
#' @return A list containing the following fields:
#' \item{cluster}{A vector of levels 1:k, indicating the cluster to which each curve is allocated.} 
#' \item{fpcaList}{A list with the fpcaObj for each separate cluster.} 
#' \item{iterToConv}{A number indicating how many iterations where required until convergence.} 
#' 
#' @examples
#' set.seed(1)
#' n <- 251
#' pts <- seq(0, 1, by=0.01)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 101) 
#' fkcObj <- FkC(sampWiener$Ly, sampWiener$Lt)
#' @references
#' \cite{Jeng-Min Chiou, Pai-Ling Li, "Functional clustering and identifying substructures of longitudinal data." Journal of the Royal Statistical Society 69 (2007): 679-699}
#' @export

FkC = function(y, t, k = 3, maxIter = 20, optns = list(maxK = 3, lean = TRUE)){ 
  
  if( (k <2) || (ceiling(length(y)*0.5) < k) ){
    warning("The value of 'k' is outside [2, 0.5*N]; reseting to 3.")
  } 
  if(is.null(optns$maxK)){
    stop("User provided 'optns' has to provided 'maxK' information.")
  }
  
  # First FPCA
  fpcaObjY <- FPCA(y, t, optns)
  
  if( fpcaObjY$optns$dataType != 'Dense' ){
    stop(paste0("The data has to be 'Dense' for FVPA to be relevant; the current dataType is : '", fpcaObjY$optns$dataType,"'!") )
  }
  
  # Initial clustering and cluster-associated FPCAs
  initialClustering <-  kmeans( fpcaObjY$xiEst, centers = k, algorithm = "Lloyd", iter.max = maxIter)
  clusterIds <- as.factor(initialClustering$cluster)
  indClustIds <- lapply(levels(clusterIds), function(u) which(clusterIds == u) )
  listOfFPCAobjs <- lapply(indClustIds, function(u) FPCA(y[u], t[u], optns) )
  
  # Iterative clustering
  oldClustIds <- clusterIds
  ymat = List2Mat(y,t); 
  convIter = FALSE
  
  for(j in 1:maxIter){ 
    newCosts <- sapply(listOfFPCAobjs, function(u) GetISEfromFPCA(u, ymat))
    newClustIds <- as.factor(apply(newCosts, 1, which.min))
   
    pendingCurves <- (sum(!(newClustIds == oldClustIds))) # Uncomment to see how many curves change cluster
    if( all(newClustIds == oldClustIds) ){
      convIter <- TRUE
      break;
    } else {
      oldClustIds <- newClustIds
      indClustIds <- lapply(levels(oldClustIds), function(u) which(oldClustIds == u) )
      listOfFPCAobjs <- lapply(indClustIds, function(u) FPCA(y[u], t[u], optns) )
    }  
  }
  
  if(!convIter){
    warning(paste0( 'FkC did not converge after maxIter = ', maxIter, ' iterations. ', pendingCurves, ' curve(s) are undecided.'))
  }
  
  return( list(cluster = newClustIds, fpcaList = listOfFPCAobjs, iterToConv = j-1) )
}  

GetISEfromFPCA = function(fpcaObj,ymat){
  numIntResults <- GetINScores(ymat , fpcaObj$obsGrid, fpcaObj$optns, fpcaObj$mu, fpcaObj$lambda, fpcaObj$phi, fpcaObj$sigma2)
  return( apply((numIntResults[['fittedY']] - ymat)^2, 1, function(y) trapzRcpp(X = fpcaObj$obsGrid, Y = y))  )
}
