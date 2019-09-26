### Functional Optimal Designs: modified cross validation for ridge parameter selection
MCVOptRidge <- function(y, t, Resp, p, RidgeCand, isDense = FALSE,
                        isRegression = FALSE, isSequential = TRUE){
  # Modified CV approach for ridge parameter selection
  n <- length(y)
  ridgecrit <- matrix(NA,ncol=2,nrow=length(RidgeCand))
  for(cviter in 1:length(RidgeCand)){ # cv loop
    ngrouping = 10
    if(isDense){ngrouping = 2}
    crittemp <- matrix(NA, ncol=2, nrow=ngrouping)
    for(ii in 1:ngrouping){
      idxA <- sort(sample(1:n, round(n/4))) # for FPCA and model component est
      idxB <- sort((1:n)[-idxA]) # for evaluation of optimal design
      yA <- y[idxA]; tA <- t[idxA]
      yB <- y[idxB]; tB <- t[idxB]
      
      obsGridA = signif(sort(unique(c(unlist(tA)))),8);
      
      nRegGridA = as.integer(1+diff(range(obsGridA))/min(diff(obsGridA))); 
      
      optnsA <- SetOptions(y = yA, t = tA, optns = list(usergrid= FALSE,nRegGrid = nRegGridA))
      CheckOptions(t = tA, optns = optnsA, n = length(yA))
      
      FPCAA <- FPCA(yA, tA, optnsA)
      CovA <- FPCAA$fittedCov; muA <- FPCAA$mu
      workGridA <- signif(FPCAA$workGrid,8);
      
      DesPoolB <- matrix(NA, ncol = 0, nrow = p)
      n_DesB <- length(idxB)
      if(isDense){n_DesB <- 1}
      for(j in 1:n_DesB){
        if(length(tB[[j]]) >= p){
          subidxj <- which(workGridA %in% tB[[j]])
          if(length(subidxj) < p){message(subidxj)}
          DesPoolB <- cbind(DesPoolB, utils::combn(subidxj,p))
        }
      }
      DesPoolB <- t(unique(t(DesPoolB))) # take only unique designs
      if(isRegression == TRUE){
        # Get CCovA
        RespA <- Resp[idxA]; RespB <- Resp[idxB]
        yA1 <- list()
        for(l in 1:length(idxA)){yA1[[l]] <- yA[[l]] * RespA[l]}
        FPCAAcross <- FPCA(yA1, tA, optnsA)
        CCovAtemp <- FPCAAcross$mu - muA*mean(RespA)
        CCovA = CCovAtemp
        #CCovA <- ConvertSupport(fromGrid = signif(FPCAAcross$obsGrid,8), toGrid = workGridA,
        #                        mu = CCovAtemp)
        # then find optdes for SR
        optdesii <- BestDes_SR_CV(p = p, ridge = RidgeCand[cviter], DesPool = DesPoolB,
                                  workGrid = workGridA, Cov = CovA,
                                  CCov = CCovA, isDense = isDense, isSequential = isSequential)
        # check out the subjects that has measurements at selected designs
        idxB2Prd <- rep(0, length(idxB))
        for(ll in 1:length(idxB)){
          if(sum(workGridA[optdesii] %in% tB[[ll]]) == p) idxB2Prd[ll] <- 1
        }
        #if(!isDense){
        #message(paste("idxB2Prd is", length(which(idxB2Prd==1)),"..."))
        #}
        yB2Prd <- yB[which(idxB2Prd == 1)]
        tB2Prd <- tB[which(idxB2Prd == 1)]
        RespB2Prd <- RespB[which(idxB2Prd == 1)]
        # make predictions
        FittedB <- RespPrd_CV(design = optdesii, ylist2Prd = yB2Prd, tlist2Prd = tB2Prd, 
                              mu = FPCAA$mu, RespTrain = RespA,
                              obsGrid = obsGridA, workGrid = workGridA,
                              ridge = RidgeCand[cviter], Cov = CovA, CCov = CCovA)
        # calculate performance criterion
        crittemp[ii,] <- APE_CV(RespB2Prd, FittedB$yres)
      } else { # for TR
        optdesii <- BestDes_TR_CV(p = p, ridge = RidgeCand[cviter], DesPool = DesPoolB,
                                  workGrid = workGridA, Cov = CovA, isDense = isDense,
                                  isSequential = isSequential)
        # check out the subjects that has measurements at selected designs
        idxB2Prd <- rep(0, length(idxB))
        for(ll in 1:length(idxB)){
          if(sum(workGridA[optdesii] %in% tB[[ll]]) == p) idxB2Prd[ll] <- 1
        }
        #if(!isDense){
        #message(paste("idxB2Prd is", length(which(idxB2Prd==1)),"..."))
        #}
        yB2Prd <- yB[which(idxB2Prd == 1)]
        tB2Prd <- tB[which(idxB2Prd == 1)]
        FittedB <- TrajRec_CV(design = optdesii, ylist2Prd = yB2Prd, tlist2Prd = tB2Prd,
                              mu = FPCAA$mu, obsGrid = obsGridA, workGrid = workGridA,
                              ridge = RidgeCand[cviter], Cov = CovA)
        # returns the fitted list at observed time points
        # calculate performance criterion
        yB2PrdAdj <- list()
        for(rr in 1:length(yB2Prd)){
          ttemp <- which( ( signif(tB2Prd[[rr]],8) %in% signif(workGridA,8) ) == TRUE)
          yB2PrdAdj[[rr]] <- yB2Prd[[rr]][ttemp]
        }
        crittemp[ii,] <- ARE_CV(yB2PrdAdj, FittedB$yres)
      }
    }
    ridgecrit[cviter,] <- colMeans(crittemp, na.rm = TRUE)
    message(paste("CV for",cviter,"out of",length(RidgeCand), "ridge parameters complete!"))
  }
  optridge = RidgeCand[which(ridgecrit[,1] == min(ridgecrit[,1]))]
  return(list(optridge = optridge, RidgeError = cbind(RidgeCand,ridgecrit)))
}

BestDes_TR_CV <- function(p, ridge, DesPool, workGrid, Cov, isDense = TRUE, isSequential = FALSE){
  # select the optimal designs for trajectory recovery case in cv
  if(isDense){
    best <- BestDes_TR(p, ridge, workGrid, Cov, isSequential=isSequential)$best
  }
  else { # sparse case
    comblist = DesPool
    temps <- rep(0,ncol(comblist))
    for(i in 1:ncol(comblist)){  temps[i] <- TRCri(comblist[,i], ridge, Cov, workGrid)  }
    best <- sort(comblist[,min(which(temps==max(temps)))])
    return(best)
    # sequential method not used for sparse case
  }
}

BestDes_SR_CV <- function(p, ridge, DesPool, workGrid, Cov, CCov, isDense = TRUE, isSequential = FALSE){
  # select the optimal designs for regression case in cv
  if(isDense){
    best <- BestDes_SR(p, ridge, workGrid, Cov, CCov, isSequential=isSequential)$best
  } else {
    comblist <- DesPool
    temps <- rep(0,ncol(comblist))
    for(i in 1:ncol(comblist)){  temps[i] <- SRCri(comblist[,i], ridge, Cov, CCov)  }
    best <- sort(comblist[,min(which(temps==max(temps)))])
    return(best)
    # sequential method not used for sparse case
  }
}

TrajRec_CV <- function(design, ylist2Prd, tlist2Prd, mu, obsGrid, workGrid, ridge, Cov){
  # recover the trajectories based on a given design vector (for cv)
  if(length(ylist2Prd) == 0){ return(list(yres = list(), tres = list()))}
  ridgeCov <- Cov + diag(ridge, nrow(Cov))
  if(length(mu) != length(workGrid)){
    mu <- ConvertSupport(fromGrid = obsGrid, toGrid = workGrid, mu = mu)
  }
  yres <- list()
  for(i in 1:length(ylist2Prd)){
    obsidx <- which((signif(workGrid,digits=8) %in% signif(tlist2Prd[[i]],8)) == TRUE)
    obsdesidx <- which( (signif(tlist2Prd[[i]],8) %in% signif(workGrid[design],8) ) == TRUE )
    yres[[i]] <- mu[obsidx] + Cov[obsidx,design] %*% solve(ridgeCov[design,design]) %*% 
      ((ylist2Prd[[i]][obsdesidx]) - mu[design])
  }
  return(list(yres = yres, tres = tlist2Prd))
}

RespPrd_CV <- function(design, ylist2Prd, tlist2Prd, mu, RespTrain, obsGrid,
                       workGrid, ridge, Cov, CCov){
  # predict the response based on a given design vector (for cv)
  if(length(ylist2Prd) == 0){ return(list(yres = c(), tres = list()))}
  ridgeCov <- Cov + diag(ridge,nrow(Cov))
  if(length(mu) != length(workGrid) || length(CCov) != length(workGrid)){
    mu <- ConvertSupport(fromGrid = obsGrid, toGrid = workGrid, mu = mu)
    CCov <- ConvertSupport(fromGrid = obsGrid, toGrid = workGrid, mu = CCov)
  }
  Ybar <- mean(RespTrain)
  res <- rep(NA, length(ylist2Prd))
  for(i in 1:length(ylist2Prd)){
    obsdesidx <- which( ( signif(tlist2Prd[[i]],8) %in% signif(workGrid[design],8)) == TRUE)
    res[i] <- Ybar + CCov[design] %*% solve(ridgeCov[design,design]) %*%
      (ylist2Prd[[i]][obsdesidx] - mu[design])
  }
  return(list(yres = res, tres = tlist2Prd))
}

ARE_CV <- function(ObsTraj, RecTraj){
  # error cirterion: MISE for trajectory recovery
  if(length(RecTraj) == 0){return(c(NA,NA))}
  nn <- length(ObsTraj)
  addi <- rep(NA, nn)
  denom <- rep(NA,nn)
  for(i in 1:nn){
    addi[i] <- sqrt(mean((ObsTraj[[i]]-RecTraj[[i]])^2))
    denom[i] <- sqrt(mean(ObsTraj[[i]]^2))
  }
  return(c(mean(addi),mean(addi)/mean(denom)))
}

APE_CV <- function(ObsResp, PredResp){
  # error criterion: MSE for regression case
  if(length(PredResp)==0){return(c(NA,NA))}
  return(c(sqrt(mean((ObsResp - PredResp)^2)),
           sqrt(mean((ObsResp - PredResp)^2))/sqrt(mean(ObsResp^2))))  
}
