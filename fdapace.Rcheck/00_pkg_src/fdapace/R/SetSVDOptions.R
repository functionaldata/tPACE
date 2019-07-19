### SetSVDOptions

SetSVDOptions <- function(Ly1, Lt1, Ly2, Lt2, SVDoptns){
  dataType1 = SVDoptns[['dataType1']]
  dataType2 = SVDoptns[['dataType2']]
  bw1 = SVDoptns[['bw1']]
  bw2 = SVDoptns[['bw2']]
  userMu1 = SVDoptns[['userMu1']]
  userMu2 = SVDoptns[['userMu2']]
  useGAM = SVDoptns[['useGAM']]
  methodSelectK = SVDoptns[['methodSelectK']]
  FVEthreshold = SVDoptns[['FVEthreshold']]
  maxK = SVDoptns[['maxK']]
  kernel = SVDoptns[['kernel']]
  nRegGrid1 = SVDoptns[['nRegGrid1']]
  nRegGrid2 = SVDoptns[['nRegGrid2']]
  bwRoutine = SVDoptns[['bwRoutine']]
  rmDiag = SVDoptns[['rmDiag']]
  noScores = SVDoptns[['noScores']]
  regulRS = SVDoptns[['regulRS']]
  flip = SVDoptns[['flip']]
  
  if(is.null(dataType1)){# do we have dataType or sparse functional data for the first sample
    dataType1 = IsRegular(Lt1);
  }
  if(is.null(dataType2)){# do we have dataType or sparse functional data for the first sample
    dataType2 = IsRegular(Lt2);
  }
  if(is.null(maxK)){ # maximum number of singular components to consider
    maxK = min(20, length(Ly1)-2)
  }
  if(is.null(kernel)){ # only gauss is avaiable for CrCov smoothing now
    kernel = "gauss"
  }
  if(is.null(nRegGrid1)){
    nRegGrid1 = 51 # currently only 51, since it is the case for GetCrCovXY()
  }
  if(is.null(nRegGrid2)){
    nRegGrid2 = 51 # currently only 51, since it is the case for GetCrCovXY()
  }
  if(is.null(userMu1)){
    userMu1 = NULL
  }
  if(is.null(userMu2)){
    userMu2 = NULL
  }
  if(is.null(bwRoutine)){
    bwRoutine = 'l-bfgs-b'
  }
  if(is.null(useGAM)){
    useGAM = FALSE
  }
    if(is.null(rmDiag)){ # whether to remove diagonal raw cov for cross cov estimation
    rmDiag = FALSE
  }
  if(is.null(FVEthreshold)){  # Default Value for the Fraction-of-Variance-Explained
    FVEthreshold = 0.9999;
  }
  if(is.null(methodSelectK)){ # the method of choosing the number of singular components K
    methodSelectK = "FVE";
  }
  if(is.numeric(methodSelectK)){
    FVEthreshold <- 1 # disable FVE selection.
    if(methodSelectK > maxK){ # check if a reasonable number of singular functions is requested
      message(paste("maxK can only be less than or equal to", maxK,"! Reset to be", maxK, "now!\n"));
    }else if(methodSelectK <= 0 || as.integer(methodSelectK)!=methodSelectK){ # check if a positive number of eigenfunctions is requested
      message("methodSelectK must be a positive integer! Reset to FVE now!\n");
      methodSelectK = "FVE"
      FVEthreshold = 0.95;
    }
  }
  
  if(is.null(noScores)){ # the method of choosing the number of singular components K
    noScores = FALSE;
  }
  if(is.null(regulRS)){ # the method of choosing the number of singular components K
    regulRS = 'sigma2';
  }
  if(is.null(flip)){ # the method of choosing the number of singular components K
    flip = FALSE;
  }
  
  
  retSVDOptns <- list(dataType1 = dataType1, dataType2 = dataType2,
                      bw1 = bw1, bw2 = bw2, userMu1 = userMu1, userMu2 = userMu2,
                      useGAM = useGAM, methodSelectK = methodSelectK,
                      FVEthreshold = FVEthreshold, maxK = maxK, flip = flip, 
                      kernel = kernel, nRegGrid1 = nRegGrid1, nRegGrid2 = nRegGrid2,
                      bwRoutine = bwRoutine, rmDiag = rmDiag, noScores = noScores, regulRS = regulRS)
  
  invalidNames <- !names(SVDoptns) %in% names(retSVDOptns)
  if (any(invalidNames)) {
    stop(sprintf('Invalid option names: %s',
                 paste0(names(SVDoptns)[invalidNames], collapse=', ')))
  }
  return(retSVDOptns)
}

