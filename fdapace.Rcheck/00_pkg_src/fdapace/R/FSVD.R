#' Functional Singular Value Decomposition
#' 
#' FSVD for a pair of dense or sparse functional data. 
#' 
#' @param Ly1 A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt1 A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param Ly2 A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt2 A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param FPCAoptns1 A list of options control parameters specified by \code{list(name=value)} for the FPC analysis of sample 1. See `?FPCA'.
#' @param FPCAoptns2 A list of options control parameters specified by \code{list(name=value)} for the FPC analysis of sample 2. See `?FPCA'.
#' @param SVDoptns A list of options control parameters specified by \code{list(name=value)} for the FSVD analysis of samples 1 & 2. See `Details`.
#'
#' @details Available control options for SVDoptns are: 
#' \describe{
#' \item{bw1}{The bandwidth value for the smoothed cross-covariance function across the direction of sample 1; positive numeric - default: determine automatically based on 'methodBwCov'}
#' \item{bw2}{The bandwidth value for the smoothed cross-covariance function across the direction of sample 2; positive numeric - default: determine automatically based on 'methodBwCov'}
#' \item{methodBwCov}{The bandwidth choice method for the smoothed covariance function; 'GMeanAndGCV' (the geometric mean of the GCV bandwidth and the minimum bandwidth),'CV','GCV' - default: 10\% of the support}
#' \item{userMu1}{The user defined mean of sample 1 used to centre it prior to the cross-covariance estimation. - default: determine automatically based by the FPCA of sample 1}
#' \item{userMu2}{The user defined mean of sample 2 used to centre it prior to the cross-covariance estimation. - default: determine automatically based by the FPCA of sample 2}
#' \item{maxK}{The maximum number of singular components to consider; default: min(20, N-1), N:# of curves.}
#' \item{kernel}{Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss"; dense data are assumed noise-less so no smoothing is performed.}
#' \item{rmDiag}{Logical describing if the routine should remove diagonal raw cov for cross cov estimation (default: FALSE) }
#' \item{noScores}{Logical describing if the routine should return functional singular scores or not (default: TRUE) }
#' \item{regulRS}{String describing if the regularisation of the compositie cross-covariance matrix should be done using 'sigma1' or 'rho' (see ?FPCA for details) (default: 'sigma2') }
#' \item{bwRoutine}{String specifying the routine used to find the optimal bandwidth 'grid-search', 'bobyqa', 'l-bfgs-b' (default: 'l-bfgs-b')}
#' \item{flip}{Logical describing if the routine should flip the sign of the singular components functions or not after the SVD of the cross-covariance matrix. (default: FALSE)}
#' \item{useGAM}{Indicator to use gam smoothing instead of local-linear smoothing (semi-parametric option) (default: FALSE)}
#' \item{dataType1}{The type of design we have for sample 1 (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV' - default:  determine automatically based on 'IsRegular'}
#' \item{dataType2}{The type of design we have for sample 2 (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV' - default:  determine automatically based on 'IsRegular'}
#' }
#' 
#' @return A list containing the following fields:
#' \item{bw1}{The selected (or user specified) bandwidth for smoothing the cross-covariance function across the support of sample 1.}
#' \item{bw2}{The selected (or user specified) bandwidth for smoothing the cross-covariance function across the support of sample 2.}
#' \item{CrCov}{The smoothed cross-covariance between samples 1 & 2.}
#' \item{sValues}{A list of length \emph{nsvd}, each entry containing the singuar value estimates for the FSC estimates.}
#' \item{nsvd}{The number of singular componentes used.}
#' \item{canCorr}{The canonical correlations for each dimension.}
#' \item{FVE}{A percentage indicating the total variance explained by chosen FSCs with corresponding 'FVEthreshold'.}
#' \item{sFun1}{An nWorkGrid by \emph{K} matrix containing the estimated singular functions for sample 1.}
#' \item{sFun2}{An nWorkGrid by \emph{K} matrix containing the estimated singular functions for sample 2.}
#' \item{grid1}{A vector of length nWorkGrid1. The internal regular grid on which the singular analysis is carried on the support of sample 1.}
#' \item{grid2}{A vector of length nWorkGrid2. The internal regular grid on which the singular analysis is carried on the support of sample 2.}
#' \item{sScores1}{A \emph{n} by \emph{K} matrix containing the singular scores for sample 1.}
#' \item{sScores2}{A \emph{n} by \emph{K} matrix containing the singular scores for sample 2.}
#' \item{optns}{A list of options used by the SVD and the FPCA's procedures.}
#' 
#' @export

FSVD <- function(Ly1, Lt1, Ly2, Lt2, FPCAoptns1 = NULL, FPCAoptns2 = NULL, SVDoptns = list()){
  # Check and refine data
  CheckData(Ly1, Lt1)
  CheckData(Ly2, Lt2)
  inputData1 <- HandleNumericsAndNAN(Ly1,Lt1);
  ly1 <- inputData1$Ly;
  lt1 <- inputData1$Lt;
  inputData2 <- HandleNumericsAndNAN(Ly2,Lt2);
  ly2 <- inputData2$Ly;
  lt2 <- inputData2$Lt;
  
  # set up and check options
  numOfCurves = length(ly1)
  SVDoptns = SetSVDOptions(Ly1=ly1, Lt1=lt1, Ly2=ly2, Lt2=lt2, SVDoptns=SVDoptns)
  #CheckSVDOptions(Ly1=y1, Lt1=t1, Ly2=y2, Lt2=t2, SVDoptns=SVDoptns)
  
  if(is.null(FPCAoptns1)){
    FPCAoptns1 = SetOptions(ly1, lt1, optns = list(dataType = SVDoptns$dataType1, nRegGrid = SVDoptns$nRegGrid1))
    CheckOptions(lt1, FPCAoptns1, numOfCurves)
  }
  
  if(is.null(FPCAoptns2)){
    FPCAoptns2 = SetOptions(ly2, lt2, optns = list(dataType = SVDoptns$dataType2, 
                                                   nRegGrid = SVDoptns$nRegGrid2))
    CheckOptions(lt2, FPCAoptns2, numOfCurves)
  }
  
  # Obtain regular grid for both samples
  obsGrid1 = sort(unique( c(unlist(lt1))))
  #regGrid1 = seq(min(obsGrid1), max(obsGrid1),length.out = SVDoptns$nRegGrid1)
  obsGrid2 = sort(unique( c(unlist(lt2))))
  #regGrid2 = seq(min(obsGrid2), max(obsGrid2),length.out = SVDoptns$nRegGrid2)
  
  ### FPCA on both functional objects
  FPCAObj1 = FPCA(Ly = ly1, Lt = lt1, optns = FPCAoptns1)
  Ymu1 = ConvertSupport(fromGrid = FPCAObj1$workGrid, toGrid = FPCAObj1$obsGrid, mu = FPCAObj1$mu)
  FPCAObj2 = FPCA(Ly = ly2, Lt = lt2, optns = FPCAoptns2)
  Ymu2 = ConvertSupport(fromGrid = FPCAObj2$workGrid, toGrid = FPCAObj2$obsGrid, mu = FPCAObj2$mu)
  if(!is.null(SVDoptns$userMu1)){
    Ymu1 = SVDoptns$userMu1
  }
  if(!is.null(SVDoptns$userMu2)){
    Ymu2 = SVDoptns$userMu2
  }
  
  ### Estimate cross covariance
  if(SVDoptns$dataType1 != 'Dense' || SVDoptns$dataType2 != 'Dense'){ # sparse or missing data
    CrCovObj = GetCrCovYX(bw1 = SVDoptns$bw1, bw2 = SVDoptns$bw2, Ly1 = ly1, Lt1 = lt1,
                          Ymu1 = Ymu1, Ly2 = ly2, Lt2 = lt2, Ymu2 = Ymu2, 
                          useGAM = SVDoptns$useGAM, rmDiag=SVDoptns$rmDiag, kern=SVDoptns$kernel,
                          bwRoutine = SVDoptns$bwRoutine)
    CrCov = CrCovObj$smoothedCC
    grid1 = CrCovObj$smoothGrid[,1]
    grid2 = CrCovObj$smoothGrid[,2]
    bw1 <- CrCovObj$bw[1]
    bw2 <- CrCovObj$bw[2]
  } else { # dense data, use sample cross covariance
    y1mat = matrix(unlist(ly1), nrow = numOfCurves, byrow = TRUE)
    y2mat = matrix(unlist(ly2), nrow = numOfCurves, byrow = TRUE)
    mudf = 0
    if(!is.null(SVDoptns$userMu1)){
      mudf = mudf + 1
    }
    if(!is.null(SVDoptns$userMu2)){
      mudf = mudf + 1
    }
    CrCov = (t(y1mat) - Ymu1) %*% t(t(y2mat) - Ymu2)/(numOfCurves - 2 + mudf)
    grid1 = obsGrid1
    grid2 = obsGrid2
    bw1 = NULL # not used for dense data
    bw2 = NULL
  }
  gridSize1 = grid1[2] - grid1[1]
  gridSize2 = grid2[2] - grid2[1]
  
  # SVD on smoothed cross covariance
  SVDObj <- svd(CrCov)
  sValues<- SVDObj$d
  sValues<- sValues[sValues>0] # pick only positive singular values
  if(length(sValues) > SVDoptns$maxK){ # reset number of singular component as maxK
    sValues<- sValues[seq_len(SVDoptns$maxK)]
  }
  # determine number of singular component retained
  if(is.numeric(SVDoptns$methodSelectK)){ # user specified number of singular components
    nsvd = SVDoptns$methodSelectK
    if(length(sValues) < nsvd){
      nsvd <- length(sValues)
      warning(sprintf("Only %d Singular Components are available, pre-specified nsvd set to %d. \n", 
                      length(sValues), length(sValues)))
    }
  } else { # select nsvd using FVEthreshold
    nsvd <- min(which(cumsum(sValues^2)/sum(sValues^2) >= SVDoptns$FVEthreshold), SVDoptns$maxK)
  }
  
  FVE <- sum(sValues[seq_len(nsvd)]^2)/sum(sValues^2) # fraction of variation explained
  sValues<- sValues[seq_len(nsvd)]
  # singular functions
  sfun1 <- SVDObj$u[,seq_len(nsvd)]
  sfun2 <- SVDObj$v[,seq_len(nsvd)]
  # normalizing singular functions
  sFun1 <- apply(sfun1, 2, function(x) {
    x <- x / sqrt(trapzRcpp(grid1, x^2))   
      return(x* ifelse(SVDoptns$flip, -1, 1))
  })
  sFun2 <- apply(sfun2, 2, function(x) {
    x <- x / sqrt(trapzRcpp(grid2, x^2))  
    return(x * ifelse(SVDoptns$flip, -1, 1))
  })
  # Grid size correction
  sValues<- sqrt(gridSize1 * gridSize2) * sValues
  
  ### calculate canonical correlations
  canCorr = rep(0, nsvd)
  for(i in seq_len(nsvd)){
    canCorr[i] = sValues[i] * (gridSize1^2 * t(sFun1[,i]) %*% FPCAObj1$fittedCov %*% (sFun1[,i]))^(-0.5) * 
      (gridSize2^2 * t(sFun2[,i]) %*% FPCAObj2$fittedCov %*% (sFun2[,i]))^(-0.5)
  }
  canCorr <- ifelse(abs(canCorr)>1,1,abs(canCorr)) * sign(canCorr) # confine this to [-1,1]
  
  if(!SVDoptns$noScores){
    
    ### calculate the singular scores
    if(SVDoptns$dataType1 != 'Dense' || SVDoptns$dataType2 != 'Dense'){ 
      # sparse data (not both dense): conditional expectation
      # patch the covariance matrix for stacked processes
      StackSmoothCov = rbind(cbind(FPCAObj1$smoothedCov, CrCov),
                             cbind(t(CrCov), FPCAObj2$smoothedCov))
      eig = eigen(StackSmoothCov)
      positiveInd <- eig[['values']] >= 0
      if (sum(positiveInd) == 0) {
        stop('All eigenvalues are negative. The covariance estimate is incorrect.')
      }
      d <- eig[['values']][positiveInd]
      eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]
      StackFittedCov <- eigenV %*% diag(x=d, nrow = length(d)) %*% t(eigenV)
      # regularize with rho/sigma2 in FPCA objects
      if( SVDoptns$regulRS == 'sigma2'){
        samp1minRS2 <- FPCAObj1$sigma2;
        samp2minRS2 <- FPCAObj2$sigma2;       
        if(is.null(samp1minRS2) ||is.null(samp2minRS2)){
          stop("Check the FPCA arguments used. At least one of the FPCA objects created does not have 'sigma2'!")
        }
      } else { #( SVDoptns$regulRS == 'rho' )
        samp1minRS2 <- FPCAObj1$rho;
        samp2minRS2 <- FPCAObj2$rho;
        if(is.null(samp1minRS2) ||is.null(samp2minRS2)){
          stop("Check the FPCA arguments used. At least one of the FPCA objects created does not have 'rho'!")
        }
      }
      
      StackRegCov <- StackFittedCov + diag(x = c(rep(samp1minRS2,length(FPCAObj1$workGrid)),
                                                 rep(samp2minRS2,length(FPCAObj2$workGrid))), 
                                           nrow = length(FPCAObj1$workGrid)+length(FPCAObj2$workGrid) )
      StackworkGrid12 = c(FPCAObj1$workGrid, FPCAObj2$workGrid + max(FPCAObj1$workGrid) + gridSize1)
      StackobsGrid12 = c(FPCAObj1$obsGrid, FPCAObj2$obsGrid + max(FPCAObj1$obsGrid) + gridSize1)
      StackRegCovObs = ConvertSupport(fromGrid = StackworkGrid12, toGrid = StackobsGrid12,
                                      Cov = StackRegCov)
      Pim = list(); Qik = list()
      Sigmai1 = list(); Sigmai2 = list()
      fittedCovObs1 = ConvertSupport(fromGrid = FPCAObj1$workGrid, toGrid = FPCAObj1$obsGrid, Cov = FPCAObj1$fittedCov)
      sfObs1 = ConvertSupport(fromGrid = grid1, toGrid = FPCAObj1$obsGrid, phi = sFun1)
      fittedCovObs2 = ConvertSupport(fromGrid = FPCAObj2$workGrid, toGrid = FPCAObj2$obsGrid, Cov = FPCAObj2$fittedCov)
      sfObs2 = ConvertSupport(fromGrid = grid2, toGrid = FPCAObj2$obsGrid, phi = sFun2)
      
      # declare singular components
      sScores1 = matrix(0, nrow = numOfCurves, ncol = nsvd)
      sScores2 = matrix(0, nrow = numOfCurves, ncol = nsvd)
      for(i in seq_len(numOfCurves)){ # calculate singular component for each obs pair
        Pim[[i]] = matrix(0, nrow=length(lt1[[i]]), ncol=nsvd)
        Qik[[i]] = matrix(0, nrow=length(lt2[[i]]), ncol=nsvd)
        Tind1 = which(FPCAObj1$obsGrid %in% lt1[[i]])
        Tind2 = which(FPCAObj2$obsGrid %in% lt2[[i]])
        Tind12 = c(Tind1, length(FPCAObj1$obsGrid) + Tind2)
        for(j in seq_len(nsvd)){
          Pim[[i]][,j] = apply(fittedCovObs1[,Tind1, drop=FALSE], 2, function(x){
            pij = trapzRcpp(X = FPCAObj1$obsGrid, Y = x*sfObs1[,j])
          })
          Qik[[i]][,j] = apply(fittedCovObs2[,Tind2, drop=FALSE], 2, function(x){
            qij = trapzRcpp(X = FPCAObj2$obsGrid, Y = x*sfObs2[,j])
          }) 
        }
        Sigmai1[[i]] = rbind(Pim[[i]], t(sValues* t(sfObs2))[Tind2, ])
        Sigmai2[[i]] = rbind(t(sValues* t(sfObs1))[Tind1, ], Qik[[i]])
        sScores1[i,] = c( t(Sigmai1[[i]]) %*% solve(StackRegCovObs[Tind12, Tind12], c(ly1[[i]] - Ymu1[Tind1], ly2[[i]] - Ymu2[Tind2])) )
        sScores2[i,] = c( t(Sigmai2[[i]]) %*% solve(StackRegCovObs[Tind12, Tind12], c(ly1[[i]] - Ymu1[Tind1], ly2[[i]] - Ymu2[Tind2])) )
      }
      
    } else { # both are dense, utilize GetINscores as it does the same job
      # sFun1/sFun2 are already in obsGrid1/obsGrid2 in dense case
      sScores1 = GetINScores(ymat = y1mat, t = lt1, optns = FPCAoptns1,
                             mu = Ymu1, lambda = sValues, phi = sFun1, sigma2 = FPCAObj1$sigma2)$xiEst    
      sScores2 = GetINScores(ymat = y2mat, t = lt2, optns = FPCAoptns2,
                             mu = Ymu2, lambda = sValues, phi = sFun2, sigma2 = FPCAObj2$sigma2)$xiEst
    }
    
  } else {
    sScores1 <- NULL
    sScores2 <- NULL
  }
  
  res <- list(bw1 = bw1, bw2 = bw2, CrCov = CrCov, sValues= sValues, nsvd = nsvd,
              canCorr = canCorr, FVE = FVE, sFun1 = sFun1, grid1 = grid1, sScores1 = sScores1,
              sFun2 = sFun2, grid2 = grid2, sScores2 = sScores2, 
              optns = list(SVDopts = SVDoptns, FPCA1opts = FPCAObj1$optns,  FPCA2opts = FPCAObj2$optns))
  class(res) <- 'FSVD'
  return(res)
}
