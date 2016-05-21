#' Functional clustering and identifying substructures of longitudinal data
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param k A scalar defining the number of clusters to define; default 3.
#' @param kSeed A scalar valid seed number to ensure replication; default: 123
#' @param maxIter A scalar defining the maximum number of iterations allowed; default 20, common for both the simple kmeans initially and the subsequent k-centres
#' @param optns A list of options control parameters specified by \code{list(name=value)} to be passed to FPCA; by default: "methodMuCovEst = 'smooth'" and "FVEthreshold = 0.90". See `Details in ?FPCA'.
#'
#' @return A list containing the following fields:
#' \item{cluster}{A vector of levels 1:k, indicating the cluster to which each curve is allocated.} 
#' \item{fpcaList}{A list with the fpcaObj for each separate cluster.} 
#' \item{iterToConv}{A number indicating how many iterations where required until convergence.} 
#' 
#' @examples
#' set.seed(1)
#' n <- 2511
#' pts <- seq(0, 1, by=0.01)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 101) 
#' kcfcObj <- kCFC(sampWiener$Ly, sampWiener$Lt)
#' @references
#' \cite{Jeng-Min Chiou, Pai-Ling Li, "Functional clustering and identifying substructures of longitudinal data." Journal of the Royal Statistical Society 69 (2007): 679-699}
#' @export

kCFC = function(y, t, k = 3, maxIter = 20, optns = list( methodMuCovEst = 'smooth', FVEthreshold = 0.90), kSeed = 123){ 
  
  if( (k <2) || (ceiling(length(y)*0.5) < k) ){
    warning("The value of 'k' is outside [2, 0.5*N]; reseting to 3.")
  } 
  if(maxIter <1){
    stop("Please allow at least 1 iteration.")
  }
  
  ## First FPCA
  fpcaObjY <- FPCA(y, t, optns)
  N <- length(y)
  if( fpcaObjY$optns$dataType != 'Dense' ){
    stop(paste0("The data has to be 'Dense' for kCFC to be relevant; the current dataType is : '", fpcaObjY$optns$dataType,"'!") )
  }
  
  ## Initial clustering and cluster-associated FPCAs
  # ## Cluster initialisation is NOT random, for each k cluster we use the medoid of the k-th dozen of points.
  # myMiniCenters <- fpcaObjY$xiEst[sapply(1:k, function(u) 
  #  which.min(rowSums(as.matrix(dist(fpcaObjY$xiEst[(1:12)+(u-1)*12,]))))) + seq(0, (k-1)*12,12),];
  
  if(!is.null(kSeed)){
    set.seed(kSeed)
  }
  
  initialClustering <- kmeans( fpcaObjY$xiEst, centers = k, algorithm = "MacQueen", iter.max = maxIter)
  clustConf0 <- as.factor(initialClustering$cluster)
  indClustIds <- lapply(levels(clustConf0), function(u) which(clustConf0 == u) )
  listOfFPCAobjs <- lapply(indClustIds, function(u) FPCA(y[u], t[u], optns) )
  
  ## Iterative clustering
  ymat <- List2Mat(y,t); 
  convInfo <- "None"
  clustConf <- list() 
  
  for(j in seq_len(maxIter)){ 
    
    # Get new costs and relevant cluster configuration
    iseCosts       <- sapply(listOfFPCAobjs, function(u) GetISEfromFPCA(u, ymat))
    clustConf[[j]] <- as.factor(apply(iseCosts, 1, which.min))
    
    # Check that clustering progressed reasonably
    if( (length(unique(clustConf[[j]])) < k) ||       # Still have k clusters
        min(summary(clustConf[[j]])) < 0.01 * N){     # Minimum cluster size is reasonable
      convInfo <- ifelse( length(unique(clustConf[[j]])) < k , "LostCluster", "TinyCluster")
      break;
    }
    # Check if algorithm converged
    if( (j>1) && any(sapply(clustConf[1:(j-1)], function(u) all(u == clustConf[[j]]))) ){
      convInfo <- "WeMadeIt!"
      break;
    } 
    
    indClustIds       <- lapply(levels(clustConf[[j]]), function(u) which(clustConf[[j]] == u) )
    listOfFPCAobjs    <- lapply(indClustIds, function(u) FPCA(y[u], t[u], optns) )
    curvesThatChanged <- ifelse(j > 1, sum(!( as.numeric(clustConf[[j]])  == as.numeric(clustConf[[j-1]] ))),
                                sum(!( as.numeric(clustConf[[j]])  == as.numeric(clustConf0))))
  } 
  
  if(convInfo == 'None'){
    warning(paste0( 'FkC did not converge after maxIter = ', maxIter, ' iterations. ', curvesThatChanged, ' curve(s) are undecided.'))
  }
  if(convInfo == 'TinyCluster'){
    warning(paste0("kCFC did not fully converge. It stopped because the smallest cluster has ",
                   "less than 1% of the samples' curves. Consider using a smaller number of clusters."))
  } 
  if(convInfo == 'LostCluster'){
    warning(paste0("kCFC did not fully converge. It stopped because it 'lost a cluster'. Consider using a smaller number of clusters."))
  }
  
  return( list(cluster = clustConf[[j]], fpcaList = listOfFPCAobjs, iterToConv = j-1, prevConf = clustConf, clustConf0 = clustConf0) )
}  

GetISEfromFPCA = function(fpcaObj,ymat){
  # First get the fitted curves for all the sample based on the mu/phi/lambda/sigma2
  # of 'fpcaObj' and then calculate their associated ISE; 'iseCost' is a n-dim vector.
  phiObs <- ConvertSupport(fpcaObj$workGrid, fpcaObj$obsGrid, phi=fpcaObj$phi)
  muObs  <- ConvertSupport(fpcaObj$workGrid, fpcaObj$obsGrid, mu=fpcaObj$mu)
  
  numIntResults <- GetINScores(ymat = ymat, t = fpcaObj$obsGrid, optns = fpcaObj$optns, mu = muObs, 
                               lambda = fpcaObj$lambda, phi = phiObs, sigma2 = fpcaObj$sigma2) 
  
  iseCost <- apply((numIntResults[['fittedY']] - ymat)^2, 1, function(y) trapzRcpp(X = fpcaObj$obsGrid, Y = y)) 
  return( iseCost )
}
