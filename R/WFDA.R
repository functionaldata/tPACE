#' Warped Functional DAta Analysis
#' 
#' Pairwise curve synchronization for functional data
#' 
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' @details WFDA uses a pairwise warping method to obtain the desired alignment (registration) of the random trajectories. The data has to be regular. The routine returns the aligned curves and the associated warping function. 
#' 
#' Available control options are 
#' \describe{
#' \item{choice}{Choice of estimating the warping functions ('weighted' or 'truncated'). If 'weighted' then weighted averages of pairwise warping functions are computed; the weighting is based on the inverse pairwise distances. If 'truncated' the pairs with the top 10% largest distances are truncated and the simple average of the remaining pairwise distances are used - default: 'truncated'}
#' \item{subsetting}{Pairwise warping functions are determined by using a subset of the whole sample; numeric (0,1] - default: 0.50}
#' \item{lambda}{Penalty parameter used for estimating pairwise warping functions; numeric - default : V*10^-4, where V is the average L2 norm of y-mean(y).}
#' \item{nknots}{Number of knots used for estimating pairwise warping functions; numeric - default: 3} 
#' \item{isPWL}{Indicator if the resulting warping functions should piece-wise linear, if FALSE 'nknots' is ignored and the resulting warping functions are simply monotonic; logical - default: TRUE} 
#' \item{seed}{Random seed for the selection of the subset of warping functions; numeric - default: 666}
#' }
#' @return A list containing the following fields: 
#' \item{lambda}{Penalty parameter used.}
#' \item{aligned}{Aligned curves evaluated at time 't'}
#' \item{h}{Warping functions for 't'} 
#' \item{hInv}{Inverse warping functions evaluated at 't'} 
#' \item{The mean cost associated with each curve}{costs} 
#' \item{timing}{The time required by time-warping.} 
#' @examples
#' N = 50;
#' eps = 0.05;
#' M = 51;
#' set.seed(123) 
#' me <- function(t) exp(-10*(t-0.5)^2);
#' T = seq(0,1,length.out = M) 
#' recondingTimesMat = matrix(nrow = N, ncol = M)
#' yMat = matrix(nrow = N, ncol = M)
#' 
#' for (i in 1:N){
#'   peak = runif(min = 0.3,max =  0.7,1)# rnorm(1, mean=0.5, sd= 0.133)
#'   recondingTimesMat[i,] = sort( c(0,    runif(min = 0, max = peak, n = round((M-3)/2)), 
#'                                   peak, runif(min = peak, max = 1, n = round((M-3)/2)), 1) )
#'   yMat[i,] = me(recondingTimesMat[i,])  + rnorm(M, mean=0.0, sd=  eps) 
#' }
#' 
#' Y = as.list(as.data.frame(t(yMat)))
#' X = rep(list(T),N)
#' 
#' sss =  WFDA(Y, X)
#' par(mfrow=c(1,2))
#' matplot(x= T, t(yMat), t='l'); grid()
#' matplot(x= T, t(sss$aligned), t='l'); grid() 
#' @references
#' \cite{Tang, R. and Mueller, H.G. (2008). "Pairwise curve synchronization for functional data." Biometrika 95, 875-889}
#' 
#' \cite{Tang, R. and Mueller, H.G. (2009) "Time-synchronized clustering of gene expression trajectories." Biostatistics 10, 32-45}
#' @export

WFDA = function(Ly, Lt, optns = list()){
  
  if(is.null(optns$isPWL)){
    optns$isPWL = TRUE
  }
  
  if(optns$isPWL){
    # Knot related checks
    if(is.null(optns$nknots)){
      optns$nknots = 3;
    } 
    if( !(optns$nknots %in% 1:7) ){
      stop("Number of knots should be between 1 and 7.")
    }
  }
  
  # Subsettig related checks
  if(is.null(optns$subsetting)){
    optns$subsetting = 0.50;
  } 
  if(findInterval(optns$subsetting - .Machine$double.eps, c(0,1)) != 1){
    stop("Number of knots should be above 0 and at most 1.")
  }
  
  # Averaging related checks
  if(is.null(optns$choice)){
    optns$choice = 'truncated';
  } 
  if( !(optns$choice %in% c('truncated','weighted') )){
    stop("The estimation of warping functions can only be done by 'truncated' or 'weighted' average.")
  }
  
  theCurrentRandomSeed <- .Random.seed #Store the current random seed to restore in the end of the function.
  #Restore the current random seed to the state it had in the beginning of the function.
  on.exit( .Random.seed <- theCurrentRandomSeed )
  
  # Seed related checks
  if(is.null(optns$seed)){
    optns$seed = 666;
  } 
  if(!(is.numeric(optns$seed))) {
    stop("The seed has to be numeric..")
  } 
   
  # Check the data validity for further analysis
  CheckData(Ly,Lt) 
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  optnsFPCA = SetOptions(Ly, Lt, list());
  
  if(optnsFPCA$dataType != 'Dense'){
    stop(paste0("The data has to be 'Dense' for WFDA to be relevant; the current dataType is : '",  optnsFPCA$dataType,"'!") )
  }
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(Ly);
  CheckOptions(Lt, optnsFPCA,numOfCurves)
  
  obsGrid = sort(unique( c(unlist(Lt))));
  
  ymat <- List2Mat(Ly, Lt)
  N = nrow(ymat)
  
  ## Super-norm normalization 
  maxAtTimeT = apply(ymat,2, function(u) max(abs(u)));
  maxAtTimeT[maxAtTimeT == 0] <- 1 ; #0/1 is still 0 so we are good.
  
  ymatNormalised <- ymat / rep(maxAtTimeT, each= N) # This is a very strong normalisation; check with Hans, I suspect it can lead to non-smooth artifacts.
  
  ## Mean function
  smcObj = GetMeanDense(ymatNormalised, obsGrid, optnsFPCA)
  
  # mu: the smoothed mean curve evaluated at times 'obsGrid'
  mu <- smcObj$mu
  
  if(is.null(optns$lambda)){
    Vy = sqrt(sum( apply(ymatNormalised,1, function(u) trapzRcpp(obsGrid, (u - mu)^2 ) ) )/(N-1))
    lambda = Vy*10^-4
  }
  
  numOfKcurves = min(round(optns$subsetting * (N-1)))
  hikMat   <- array(dim = c(numOfKcurves,M,N) ) 
  distMat <- matrix( nrow = N, ncol = numOfKcurves)
  hMat   <- array(dim = c(N,M) )
  hInvMat   <- array(dim = c(N,M) )
  alignedMat   <- array(dim = c(N,M) )
  
  getcurveJ <- function(tj, curvej){
    approx(x = obsGrid, y = curvej, xout = tj)$y
  }
  
  theCost <- function(curvei, curvek,lambda,ti,tk){ 
    sum((getcurveJ(tk, curvek)-curvei)^2) + lambda * sum((tk-ti)^2) 
  }
  
  getHikRS <- function(curvei, curvek, lambda){
    myCosts <- sapply(1:numOfKcurves^2, function(u){ set.seed(u); 
      theCost(curvei, curvek, lambda, obsGrid,  c(0,Rcppsort(runif(M-2)) ,1)) })
    set.seed( which.min( myCosts ) )
    minCost <- min(myCosts)  
    return( c(0,Rcppsort(runif(M-2)),1) )
  }
  
  getSol <- function(x){
    approx(x = seq(0,1, length.out = (2+ optns$nknots)), y = c(0, Rcppsort(x),1) ,n = M)$y
  }
  
  theCostOptim <- function(x , curvei, curvek,lambda,ti){
    tk = getSol(x);
    sum((getcurveJ(tk, curvek)-curvei)^2) + lambda * sum((tk-ti)^2) 
  }
  
  getHikOptim <- function(curvei, curvek, lambda, minqaAvail ){
    s0 <- seq(0,1,length.out = (2+ optns$nknots))[2:(1+optns$nknots)]
    if( !minqaAvail ) { 
      optimRes <- optim( par = s0, fn = theCostOptim, method = 'L-BFGS-B', 
                         lower = rep(1e-6, optns$nknots), upper = rep(1 - 1e-6, optns$nknots),
                         curvei = curvei, curvek = curvek, lambda = lambda, ti =obsGrid)
    } else {
      optimRes <-  minqa::bobyqa( par = s0, fn = theCostOptim,  
                                  lower = rep(1e-6, optns$nknots), upper = rep(1 - 1e-6, optns$nknots),
                                  curvei = curvei, curvek = curvek, lambda = lambda, ti =obsGrid)
    }
    bestSol <- getSol(optimRes$par)  
    return( bestSol )
  }
  
  
  start <- Sys.time ()
  if( !is.element('minqa', installed.packages()[,1]) && optns$isPWL){
    warning("Cannot use 'minqa::bobyqa' to find the optimal knot locations as 'minqa' is not installed. We will do an 'L-BFGS-B' search.") 
    minqaAvail = FALSE
  } else {
    minqaAvail = TRUE
  }
  
  for(i in seq_len(N)){ # For each curve
    # print(i)
    set.seed( i + optns$seed );
    curvei = ymatNormalised[i,];
    candidateKcurves = sample(seq_len(N)[-i], numOfKcurves)  
    
    for(k in seq_len(numOfKcurves)){ # For each of the candidate curves 
      if(!optns$isPWL){
        hikMat[k, ,i] = getHikRS(curvei, ymatNormalised[candidateKcurves[k],], lambda)
      } else {
        hikMat[k, ,i] = getHikOptim(curvei, ymatNormalised[candidateKcurves[k],], lambda, minqaAvail)
      }
      distMat[i,k] =  mean( ( getcurveJ(tj =hikMat[k, ,i], curvei) - ymatNormalised[candidateKcurves[k]])^2 )
    }
    
    if(optns$choice == 'weighted'){
      hMat[i,] = apply(hikMat[, ,i] , 2, weighted.mean, distMat[i,])
    } else {
      hMat[i,] = apply(hikMat[  (distMat[i,] <= quantile( distMat[i,], p=0.90) ),  ,i] , 2, mean)
    }
    
    hInvMat[i,] = approx(y = obsGrid, x = hMat[i,], xout = obsGrid)$y
    alignedMat[i,] = approx(x = obsGrid, y = ymat[i,], xout = hInvMat[i,])$y
    
  }  
  
  timing = Sys.time () - start
  ret <- list(lambda = lambda, h = hMat, hInv = hInvMat, aligned = alignedMat, costs = rowMeans(distMat), timing = timing)
  class(ret) <- 'WFDA' 
  
  return(ret); 
}
