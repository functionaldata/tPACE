#' Optimal Designs for Functional and Longitudinal Data
#' for Trajectory Recovery or Scalar Response Prediction
#'
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param Resp A vector of response values 
#' @param p A prefixed positive integer indicating the number of optimal design points requested, with default: 3.
#' @param optns A list of options control parameters specified by \code{list(name=value)} for FPCA, with default: list().
#' @param isRegression A logical argument, indicating the purpose of the optimal designs: TRUE for scalar response prediction, FALSE for trajectory recovery (default).
#' @param isSequential A logical argument, indicating whether to use the sequential optimization procedure for faster computation, recommended for relatively large p (default: FALSE).
#' @param RidgeCand A vector of positive ridge penalty candidates for regularization.
#' 
#' @details To select a proper RidgeCand, check with the returned optimal ridge parameter. The selected parameter is valid if not on the boundary of the range.
#' 
#' @return A list containing the following fields:
#' \item{OptDes}{The vector of optimal design points of the regular time grid of the observed data.}
#' \item{R2}{Coefficient of determination. (Check the paper for details.)}
#' \item{R2adj}{Adjusted coefficient of determination.}
#' \item{OptRidge}{The selected ridge parameter.}
#' 
#' @examples
#' set.seed(1)
#' n <- 50
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 21)
#' res <- Foptdes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, Resp=NULL, p=3, optns=list(),
#'                isSequential=TRUE, RidgeCand = seq(0.1,1,0.1))
#' @references
#' \cite{Ji, H., Mueller, H.G. (2016) "Optimal Designs for Longitudinal and Functional Data" Journal of the Royal Statistical Society: Series B (Statistical Methodology)}
#' 
#' @export

Foptdes <- function(Ly = NULL, Lt = NULL, Resp = NULL, p = 3, optns = list(),
                        isRegression = FALSE, isSequential = FALSE, RidgeCand = NULL){
  CheckData(y = Ly, t = Lt);
  inputData <- HandleNumericsAndNAN(Ly, Lt);
  y <- inputData$Ly;
  t <- inputData$Lt;  
  
  optns = SetOptions(y, t, optns);
  
  numOfCurves = length(y);
  CheckOptions(t, optns, numOfCurves);
  if(optns$dataType == "Dense"){
    isDense = TRUE;
  } else {
    isDense = FALSE; # currently dense with missing is treated as sparse
  }
  
  obsGrid = sort(unique(c(unlist(t))));
  RegGrid = seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid);

  # find the best ridge parameter via cross validation
  OptRidge <- MCVoptRidge(y = y, t = t, Resp = Resp, p = p, RidgeCand = RidgeCand,
                          isDense = isDense, nRegGrid = optns$nRegGrid,
                          isRegression = isRegression, isSequential = isSequential)
  optridge <- OptRidge$optridge
  # find optdes with optridge
  TrainFPCA <- FPCA(Ly=y, Lt=t, optns=optns)
  if(isRegression == FALSE){ # Trajectory Recovery
    BestDesTR <- bestdes.TR(p=p, ridge=optridge, workGrid=TrainFPCA$workGrid,
                            Cov=TrainFPCA$fittedCov, isSequential=isSequential)$best
    # calculate R2_X
    VarX <- sum(TrainFPCA$lambda)
    mu <- TrainFPCA$mu
    Cov <- TrainFPCA$fittedCov
    ridgeCov <- TrainFPCA$fittedCov + diag(optridge, nrow(Cov))
    R2XNum <- sum(diag(Cov[,BestDesTR] %*% solve(ridgeCov[BestDesTR, BestDesTR]) %*% Cov[BestDesTR,]))*diff(RegGrid)[1]
    R2X <- R2XNum/VarX
    R2Xadj <- 1-(1-R2X)*(length(y)-1)/(length(y)-p-1)
    return(list(OptDes = RegGrid[BestDesTR], R2 = R2X, R2adj = R2Xadj, OptRidge = OptRidge))
  } else{ # scalar response regression
    mu <- TrainFPCA$mu
    Cov <- TrainFPCA$fittedCov
    ridgeCov <- TrainFPCA$fittedCov + diag(optridge, nrow(Cov))
    # FPCA for cross cov
    y1 <- list()
    for(subj in 1:length(y)){y1[[subj]] = y[[subj]]*Resp[subj]}
    FPCACC <- FPCA(y1, t, optns)
    CCovtemp <- FPCACC$mu - mean(Resp)*mu
    CCov <- ConvertSupport(fromGrid = TrainFPCA$obsGrid, toGrid = TrainFPCA$workGrid, mu = CCovtemp)
    BestDesSR <- bestdes.SR(p=p, ridge=optridge, workGrid=TrainFPCA$workGrid,
                            Cov=TrainFPCA$fittedCov, CCov=CCov, isSequential=isSequential)$best
    R2Y <- (var(Resp) - CCov[BestDesSR] %*% solve(ridgeCov[BestDesSR,BestDesSR]) %*% CCov[BestDesSR])/var(Resp)
    R2Yadj <- 1-(1-R2Y)*(length(y)-1)/(length(y)-p-1)
    return(list(OptDes = RegGrid[BestDesSR], R2 = R2Y, R2adj = R2Yadj, OptRidge = OptRidge))
  }
}

bestdes.TR <- function(p, ridge, workGrid, Cov, isSequential=FALSE, isMed = FALSE){
  # select optimal designs for trajectory recovery case, sequential method available
  if(isSequential == FALSE){ # Global Selection
    comblist <- combn(1:length(workGrid),p)
    temps <- rep(0,ncol(comblist))
    for(i in 1:ncol(comblist)){  temps[i] <- TRcri(comblist[,i], ridge, Cov, workGrid)  }
    best <- sort(comblist[,min(which(temps==max(temps)))])
    if(isMed == TRUE) {
      med <- sort(comblist[,min(which(temps==sort(temps)[round(ncol(comblist)/2)]))])
      return(list(best=best,med=med))
    } else {return(list(best=best,med=NULL))}
  } else { # Sequential optimization
    optdes <- c()
    for(iter in 1:p){
      candidx <- which(!((1:length(workGrid)) %in% optdes))
      seqcri <- rep(NA, length(candidx))
      for(i in 1:length(candidx)){
        tempdes <- sort(c(optdes,candidx[i]))
        seqcri[i] <- TRcri(tempdes, ridge, Cov, workGrid)
      }
      optdes <- sort(c(optdes, candidx[min(which(seqcri == max(seqcri)))]))
    }
    return(list(best=optdes,med=NULL)) # based on sequential selection
  }
}  

bestdes.SR <- function(p, ridge, workGrid, Cov, CCov, isSequential=FALSE, isMed = FALSE){
  # select optimal designs for regression case, sequential method available
  if(isSequential == FALSE){
    comblist <- combn(1:length(workGrid), p)
    temps <- rep(0,ncol(comblist))
    for(i in 1:ncol(comblist)){  temps[i] <- SRcri(comblist[,i], ridge, Cov, CCov)  }
    best <- sort(comblist[,min(which(temps==max(temps)))])
    if(isMed == TRUE){
      med <- sort(comblist[,min(which(temps==sort(temps)[round(ncol(comblist)/2)]))])
      return(list(best=best,med=med))
    } else { return(list(best=best,med=NULL)) }
  } else{ # sequential selection
    optdes <- c()
    for(iter in 1:p){
      candidx <- which(!((1:length(workGrid)) %in% optdes))
      seqcri <- rep(NA, length(candidx))
      for(i in 1:length(candidx)){
        tempdes <- sort(c(optdes,candidx[i]))
        seqcri[i] <- SRcri(tempdes, ridge, Cov, CCov)
      }
      optdes <- sort(c(optdes, candidx[min(which(seqcri == max(seqcri)))]))
    }
    return(list(best=optdes,med=NULL))
  }
}

TRcri <- function(design, ridge, Cov, workGrid){
  # Optimization criterion for TR
  # Numerical Integration, equal to matrix multiplication if time grid is year
  design <- sort(design)
  RidgeCov <- Cov + diag(ridge, nrow(Cov))
  designcovinv <- solve(RidgeCov[design,design])
  if(length(design) > 1){
    trcri <- trapz(workGrid,diag(t(Cov[design,])%*%designcovinv%*%(Cov[design,])))
  } else {
    trcri <- trapz(workGrid,diag(Cov[design,]%*%designcovinv%*%(Cov[design,])))
  }
  return(trcri)
}

SRcri <- function(design,ridge,Cov,CCov){
  # Optimization criterion for SR
  design <- sort(design)
  ridgeCov <- Cov + diag(ridge,nrow(Cov))
  srcri <- t(CCov[design]) %*% solve(ridgeCov[design,design]) %*% CCov[design]
  return(srcri)
}

bestdes.TR.cv <- function(p, ridge, DesPool, workGrid, Cov, isDense = TRUE, isSequential = FALSE){
  # select the optimal designs for trajectory recovery case in cv
  if(isDense){
    best <- bestdes.TR(p, ridge, workGrid, Cov, isSequential=isSequential)$best
  }
  else { # sparse case
  comblist = DesPool
  temps <- rep(0,ncol(comblist))
  for(i in 1:ncol(comblist)){  temps[i] <- TRcri(comblist[,i], ridge, Cov, workGrid)  }
  best <- sort(comblist[,min(which(temps==max(temps)))])
  return(best)
  # sequential method not used for sparse case
  }
}

bestdes.SR.cv <- function(p, ridge, DesPool, workGrid, Cov, CCov, isSequential = FALSE){
  # select the optimal designs for regression case in cv
  if(isDense){
    best <- bestdes.SR(p, ridge, workGrid, Cov, CCov, isSequential=isSequential)$best
  } else {
  comblist <- DesPool
  temps <- rep(0,ncol(comblist))
  for(i in 1:ncol(comblist)){  temps[i] <- SRcri(comblist[,i], ridge, Cov, CCov)  }
  best <- sort(comblist[,min(which(temps==max(temps)))])
  return(best)
  # sequential method not used for sparse case
  }
}

TrajRec.cv <- function(design, ylist2Prd, tlist2Prd, mu, obsGrid, workGrid, ridge, Cov){
  # recover the trajectories based on a given design vector (for cv)
  if(length(ylist2Prd) == 0){ return(list(yres = list(), tres = list()))}
  ridgeCov <- Cov + diag(ridge, nrow(Cov))
  if(length(mu) != length(workGrid)){
    mu <- ConvertSupport(fromGrid = obsGrid, toGrid = workGrid, mu = mu)
  }
  yres <- list()
  for(i in 1:length(ylist2Prd)){
    obsidx <- which((signif(workGrid,digits=8) %in% signif(tlist2Prd[[i]],8)) == TRUE)
    obsdesidx <- which((tlist2Prd[[i]] %in% workGrid[design]) == TRUE)
    yres[[i]] <- mu[obsidx] + Cov[obsidx,design] %*% solve(ridgeCov[design,design]) %*% 
      ((ylist2Prd[[i]][obsdesidx]) - mu[design])
  }
  return(list(yres = yres, tres = tlist2Prd))
}

RespPrd.cv <- function(design, ylist2Prd, tlist2Prd, mu, RespTrain, obsGrid,
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
    obsdesidx <- which((signif(tlist2Prd[[i]],8) %in% signif(workGrid[design],8)) == TRUE)
    res[i] <- Ybar + CCov[design] %*% solve(ridgeCov[design,design]) %*%
      (ylist2Prd[[i]][obsdesidx] - mu[design])
  }
  return(list(yres = res, tres = tlist2Prd))
}

ARE.cv <- function(ObsTraj, RecTraj){
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

APE.cv <- function(ObsResp, PredResp){
  # error criterion: MSE for regression case
  if(length(PredResp)==0){return(c(NA,NA))}
  return(c(sqrt(mean((ObsResp - PredResp)^2)),
           sqrt(mean((ObsResp - PredResp)^2))/sqrt(mean(ObsResp^2))))  
}

MCVoptRidge <- function(y, t, Resp, p, RidgeCand, isDense = FALSE, nRegGrid,
                        isRegression = FALSE, isSequential = TRUE){
  # Modified CV approach for ridge parameter selection
  n <- length(y)
  ridgecrit <- matrix(NA,ncol=2,nrow=length(RidgeCand))
  for(cviter in 1:length(RidgeCand)){ # cv loop
    ngrouping = 10
    if(isDense){ngrouping = 2}
    crittemp <- matrix(NA, ncol=2, nrow=ngrouping)
    for(ii in 1:ngrouping){
      idxA <- sort(sample(1:n, round(n/2))) # for FPCA and model component est
      idxB <- sort((1:n)[-idxA]) # for evaluation of optimal design
      yA <- y[idxA]; tA <- t[idxA]
      yB <- y[idxB]; tB <- t[idxB]
      optnsA <- SetOptions(y = yA, t = tA, optns = list(nRegGrid = nRegGrid))
      CheckOptions(t = tA, optns = optnsA, n = length(yA))

      FPCAA <- FPCA(yA, tA, optnsA)
      CovA <- FPCAA$fittedCov; muA <- FPCAA$mu
      workGridA <- signif(FPCAA$workGrid,8); obsGridA <- signif(FPCAA$obsGrid,8)
      
      DesPoolB <- matrix(NA, ncol = 0, nrow = p)
      n_DesB <- length(idxB)
      if(isDense){n_DesB <- 1}
      for(j in 1:n_DesB){
        if(length(tB[[j]]) >= p){
          subidxj <- which(tB[[j]] %in% workGridA)
          DesPoolB <- cbind(DesPoolB, combn(subidxj,p))
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
        CCovA <- ConvertSupport(fromGrid = signif(FPCAAcross$obsGrid,8), toGrid = workGridA,
                                mu = CCovAtemp)
        # then find optdes for SR
        optdesii <- bestdes.SR.cv(p = p, ridge = RidgeCand[cviter], DesPool = DesPoolB,
                                    workGrid = workGridA, Cov = CovA,
                                    CCov = CCovA, isSequential = isSequential)
        # check out the subjects that has measurements at selected designs
        idxB2Prd <- rep(0, length(idxB))
        for(ll in 1:length(idxB)){
          if(sum(workGridA[optdesii] %in% tB[[ll]]) == p) idxB2Prd[ll] <- 1
        }
        #if(isDense){
        cat(paste("idxB2Prd is", length(which(idxB2Prd==1)),"...\n"))
        #}
        yB2Prd <- yB[which(idxB2Prd == 1)]
        tB2Prd <- tB[which(idxB2Prd == 1)]
        RespB2Prd <- RespB[which(idxB2Prd == 1)]
        # make predictions
        FittedB <- RespPrd.cv(design = optdesii, ylist2Prd = yB2Prd, tlist2Prd = tB2Prd, 
                                   mu = FPCAA$mu, RespTrain = RespA,
                                   obsGrid = obsGridA, workGrid = workGridA,
                                   ridge = RidgeCand[cviter], Cov = CovA, CCov = CCovA)
        # calculate performance criterion
        crittemp[ii,] <- APE.cv(RespB2Prd, FittedB$yres)
      } else { # for TR
        optdesii <- bestdes.TR.cv(p = p, ridge = RidgeCand[cviter], DesPool = DesPoolB,
                               workGrid = workGridA, Cov = CovA, isDense = isDense,
                               isSequential = isSequential)
        # check out the subjects that has measurements at selected designs
        idxB2Prd <- rep(0, length(idxB))
        for(ll in 1:length(idxB)){
          if(sum(workGridA[optdesii] %in% tB[[ll]]) == p) idxB2Prd[ll] <- 1
        }
        cat(paste("idxB2Prd is", length(which(idxB2Prd==1)),"...\n"))
        yB2Prd <- yB[which(idxB2Prd == 1)]
        tB2Prd <- tB[which(idxB2Prd == 1)]
        FittedB <- TrajRec.cv(design = optdesii, ylist2Prd = yB2Prd, tlist2Prd = tB2Prd,
                                mu = FPCAA$mu, obsGrid = obsGridA, workGrid = workGridA,
                                ridge = RidgeCand[cviter], Cov = CovA)
        # returns the fitted list at observed time points
        # calculate performance criterion
        crittemp[ii,] <- ARE.cv(yB2Prd, FittedB$yres)
      }
    }
    ridgecrit[cviter,] <- colMeans(crittemp, na.rm = TRUE)
    cat(paste("CV for",cviter,"out of",length(RidgeCand), "ridge parameters complete!\n"))
  }
  optridge = RidgeCand[which(ridgecrit[,1] == min(ridgecrit[,1]))]
  return(list(optridge = optridge, RidgeError = cbind(RidgeCand,ridgecrit)))
}
