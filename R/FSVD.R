# Functional Singular Component Analysis of two functions
# Procedures
# 1) Apply the functional singular value decomposition to the kernel
#    estimator of cross-covariance surface between X and Y.
#
# 2) Computation of the singular values and the singular functions.
#
# 3) Computation of the functional correlation. # wrapper or included with options

FSVD <- function(Ly1, Lt1, Ly2, Lt2, #FPCAoptns1 = list(), FPCAoptns2 = list(),
                 SVDoptns = list()){
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
  FPCAoptns1 = SetOptions(ly1, lt1, optns = list(dataType = SVDoptns$dataType1, 
                                                 nRegGrid = SVDoptns$nRegGrid1))
  CheckOptions(lt1, FPCAoptns1, numOfCurves)
  FPCAoptns2 = SetOptions(ly2, lt2, optns = list(dataType = SVDoptns$dataType2, 
                                                 nRegGrid = SVDoptns$nRegGrid2))
  CheckOptions(lt2, FPCAoptns2, numOfCurves)
  
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
                          quadApprox = SVDoptns$quadApprox)
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
  sv <- SVDObj$d
  sv <- sv[sv>0] # pick only positive singular values
  if(length(sv) > SVDoptns$maxK){ # reset number of singular component as maxK
    sv <- sv[1:SVDoptns$maxK]
  }
  # determine number of singular component retained
  if(is.numeric(SVDoptns$methodSelectK)){ # user specified number of singular components
    nsvd = SVDoptns$methodSelectK
    if(length(sv) < nsvd){
      nsvd <- length(sv)
      warning(sprintf("Only %d Singular Components are available, pre-specified nsvd set to %d. \n", 
                      length(sv), length(sv)))
    }
  } else { # select nsvd using FVEthreshold
    nsvd <- min(which(cumsum(sv^2)/sum(sv^2) >= SVDoptns$FVEthreshold), SVDoptns$maxK)
  }
  
  FVE <- sum(sv[1:nsvd]^2)/sum(sv^2) # fraction of variation explained
  sv <- sv[1:nsvd]
  # singular functions
  sfun1 <- SVDObj$u[,1:nsvd]
  sfun2 <- SVDObj$v[,1:nsvd]
  # normalizing singular functions
  sf1 <- apply(sfun1, 2, function(x) {
    x <- x / sqrt(trapzRcpp(grid1, x^2)) 
      return(x)
  })
  sf2 <- apply(sfun2, 2, function(x) {
    x <- x / sqrt(trapzRcpp(grid2, x^2)) 
    return(x)
  })
  # Grid size correction
  sv <- sqrt(gridSize1 * gridSize2) * sv
  
  ### calculate canonical correlations
  rho = rep(0, nsvd)
  for(i in 1:nsvd){
    rho[i] = sv[i] * (gridSize1^2 * t(sf1[,i]) %*% FPCAObj1$fittedCov %*% (sf1[,i]))^(-0.5) * 
      (gridSize2^2 * t(sf2[,i]) %*% FPCAObj2$fittedCov %*% (sf2[,i]))^(-0.5)
  }
  
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
    StackRegCov <- StackFittedCov + diag(x = c(rep(FPCAObj1$sigma2,length(FPCAObj1$workGrid)),
                                               rep(FPCAObj2$sigma2,length(FPCAObj2$workGrid))), 
                                         nrow = length(FPCAObj1$workGrid)+length(FPCAObj2$workGrid) )
    StackworkGrid12 = c(FPCAObj1$workGrid, FPCAObj2$workGrid + max(FPCAObj1$workGrid) + gridSize1)
    StackobsGrid12 = c(FPCAObj1$obsGrid, FPCAObj2$obsGrid + max(FPCAObj1$obsGrid) + gridSize1)
    StackRegCovObs = ConvertSupport(fromGrid = StackworkGrid12, toGrid = StackobsGrid12,
                                    Cov = StackRegCov)
    Pim = list(); Qik = list()
    Sigmai1 = list(); Sigmai2 = list()
    fittedCovObs1 = ConvertSupport(fromGrid = FPCAObj1$workGrid, toGrid = FPCAObj1$obsGrid, Cov = FPCAObj1$fittedCov)
    sfObs1 = ConvertSupport(fromGrid = grid1, toGrid = FPCAObj1$obsGrid, phi = sf1)
    fittedCovObs2 = ConvertSupport(fromGrid = FPCAObj2$workGrid, toGrid = FPCAObj2$obsGrid, Cov = FPCAObj2$fittedCov)
    sfObs2 = ConvertSupport(fromGrid = grid2, toGrid = FPCAObj2$obsGrid, phi = sf2)
    
    # declare singular components
    sc1 = matrix(0, nrow = numOfCurves, ncol = nsvd)
    sc2 = matrix(0, nrow = numOfCurves, ncol = nsvd)
    for(i in 1:numOfCurves){ # calculate singular component for each obs pair
      Pim[[i]] = matrix(0, nrow=length(lt1[[i]]), ncol=nsvd)
      Qik[[i]] = matrix(0, nrow=length(lt2[[i]]), ncol=nsvd)
      Tind1 = which(FPCAObj1$obsGrid %in% lt1[[i]])
      Tind2 = which(FPCAObj2$obsGrid %in% lt2[[i]])
      Tind12 = c(Tind1, length(FPCAObj1$obsGrid) + Tind2)
      for(j in 1:nsvd){
        Pim[[i]][,j] = apply(fittedCovObs1[,Tind1], 2, function(x){
          pij = trapzRcpp(X = FPCAObj1$obsGrid, Y = x*sfObs1[,j])
        })
        Qik[[i]][,j] = apply(fittedCovObs2[,Tind2], 2, function(x){
          qij = trapzRcpp(X = FPCAObj2$obsGrid, Y = x*sfObs2[,j])
        }) 
      }
      Sigmai1[[i]] = rbind(Pim[[i]], t(sv * t(sfObs2))[Tind2, ])
      Sigmai2[[i]] = rbind(t(sv * t(sfObs1))[Tind1, ], Qik[[i]])
      sc1[i,] = c( t(Sigmai1[[i]]) %*% solve(StackRegCovObs[Tind12, Tind12]) %*% 
        c(ly1[[i]] - Ymu1[Tind1], ly2[[i]] - Ymu2[Tind2]) )
      sc2[i,] = c( t(Sigmai2[[i]]) %*% solve(StackRegCovObs[Tind12, Tind12]) %*% 
        c(ly1[[i]] - Ymu1[Tind1], ly2[[i]] - Ymu2[Tind2]) )     
    }
    
  } else { # both are dense, utilize GetINscores as it does the same job
    # sf1/sf2 are already in obsGrid1/obsGrid2 in dense case
    sc1 = GetINScores(ymat = y1mat, t = lt1, optns = FPCAoptns1,
                      mu = Ymu1, lambda = sv, phi = sf1, sigma2 = FPCAObj1$sigma2)$xiEst    
    sc2 = GetINScores(ymat = y2mat, t = lt2, optns = FPCAoptns2,
                      mu = Ymu2, lambda = sv, phi = sf2, sigma2 = FPCAObj2$sigma2)$xiEst
  }
  
  res <- list(bw1 = bw1, bw2 = bw2, CrCov = CrCov, sv = sv, nsvd = nsvd,
              rho = rho, FVE = FVE, sf1 = sf1, grid1 = grid1, sc1 = sc1,
              sf2 = sf2, grid2 = grid2, sc2 = sc2, optns = SVDoptns)
  class(res) <- 'FSVD'
  return(res)
}
