#' Functional clustering and identifying substructures of longitudinal data
#' 
#' Default:  Cluster functional data using the functional principal component (FPC) scores obtained from the data
#' FPC analysis using EMCluster (Chen and Maitra, 2015) or directly clustering the functional data using kCFC (Chiou and Li, 2007).
#' 
#' Within EMCluster, uses the model initiated \code{EMCluster::em.EM} and returns the optimal model based on \code{EMCluster::emcluster}. 
#' See ?EMCluster::emcluster for details.
#' 
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param k A scalar defining the number of clusters to define; default 3.
#' @param cmethod A string specifying the clustering method to use ('EMCluster' or 'kCFC'); default: 'EMCluster'. 
#' @param optnsFPCA A list of options control parameters specified by \code{list(name=value)} to be used for by FPCA on the sample y; by default: 
#' "list( methodMuCovEst ='smooth', FVEthreshold= 0.90, methodBwCov = 'GCV', methodBwMu = 'GCV' )". See `Details in ?FPCA'.
#' @param optnsCS A list of options control parameters specified by \code{list(name=value)} to be used for cluster-specific FPCA from kCFC; by default:  
#' "list( methodMuCovEst ='smooth', FVEthreshold= 0.70, methodBwCov = 'GCV', methodBwMu = 'GCV' )". See `Details in ?FPCA' and '?kCFC'. 
#'
#' @return A list containing the following fields:
#' \item{cluster}{A vector of levels 1:k, indicating the cluster to which each curve is allocated.} 
#' \item{fpca}{An FPCA object derived from the sample used by Rmixmod, otherwise NULL.} 
#' \item{clusterObj}{Either a EMCluster object or kCFC object.} 
#' 
#' @examples
#' \donttest{
#' data(medfly25)
#' Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs) 
#' newClust <- FClust(Flies$Ly, Flies$Lt, k = 2, optnsFPCA = 
#'                     list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90))
#'                     
#' # We denote as 'veryLowCount' the group of flies that lay less
#' # than twenty-five eggs during the 25-day period examined.
#' 
#' veryLowCount = ifelse( sapply( unique(medfly25$ID), function(u) 
#'                    sum( medfly25$nEggs[medfly25$ID == u] )) < 25, 0, 1)
#' N <- length(unique(medfly25$ID))
#' (correctRate <- sum( (1 + veryLowCount) ==  newClust$cluster) / N) # 99.6%
#' }
#' @references
#' \cite{Wei-Chen Chen and Ranjan Maitra, "EMCluster: EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian Distribution". (2015)}
#'
#' \cite{Julien Jacques and Cristian Preda, "Funclust: A curves clustering method using functional random variables density approximation". Neurocomputing 112 (2013): 164-171}
#' 
#' \cite{Jeng-Min Chiou and Pai-Ling Li, "Functional clustering and identifying substructures of longitudinal data". Journal of the Royal Statistical Society B 69 (2007): 679-699}
#' @export

FClust = function(Ly, Lt, k = 3, cmethod = 'EMCluster', optnsFPCA = NULL, optnsCS = NULL){ 
  
  if(is.null(optnsFPCA)){
     optnsFPCA = list( methodMuCovEst = 'smooth', FVEthreshold = 0.90, methodBwCov = 'GCV', methodBwMu = 'GCV')
  }
  if(is.null(optnsCS)){
    optnsCS = list( methodMuCovEst = 'smooth', FVEthreshold = 0.70, methodBwCov = 'GCV', methodBwMu = 'GCV')
  }
  
  if( (k <2) || (floor(length(Ly)*0.5) < k) ){
    warning("The value of 'k' is outside [2, 0.5*N]; clustering is possibly incoherent.")
  } 
  if( !(cmethod %in% c("EMCluster", "kCFC")) ){
    stop("The clustering method specified in neither 'EMCluster' or 'kCFC'.")
  }
  
  if( cmethod == 'EMCluster'){
    if( !requireNamespace("EMCluster", quietly=TRUE) ) {#!is.element('EMCluster', installed.packages()[,1])
      stop("Cannot the use the EMCluster method; the package 'EMCluster' is unavailable.")
    }
    # suppressMessages(library(EMCluster))
    fpcaObjY <- FPCA(Ly = Ly, Lt = Lt, optnsFPCA)
    xiData <- as.data.frame(fpcaObjY$xiEst) 
    clusterObj <- EMCluster::emcluster(xiData, EMCluster::em.EM(xiData, nclass =k, EMC = EMCluster::.EMControl()), assign.class= TRUE, EMC = EMCluster::.EMControl() )
    clustConf = clusterObj$class 
  } else {
    fpcaObjY <- NULL
    clusterObj <- kCFC(y= Ly, t= Lt, k = k, optnsSW = optnsFPCA, optnsCS = optnsCS)
    clustConf <- clusterObj$cluster
  }
  
  return( list(cluster = clustConf, fpca = fpcaObjY, clusterObj =  clusterObj) )
}  
