#' Optimal Designs for Functional and Longitudinal Data
#' for Trajectory Recovery or Scalar Response Prediction
#'
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param Resp A vector of response values, keep void for trajectory recovery, only necessary for scalar response prediction task.
#' @param p A fixed positive integer indicating the number of optimal design points requested, with default: 3.
#' @param optns A list of options control parameters specified by \code{list(name=value)} for FPCA, with default: list().
#' @param isRegression A logical argument, indicating the purpose of the optimal designs: TRUE for scalar response prediction, FALSE for trajectory recovery, with default value !missing(Resp).
#' @param isSequential A logical argument, indicating whether to use the sequential optimization procedure for faster computation, recommended for relatively large p (default: FALSE).
#' @param RidgeCand A vector of positive numbers as ridge penalty candidates for regularization. The final value is selected via cross validation. If only 1 ridge parameter is specified, CV procedure is skipped.
#' 
#' @details To select a proper RidgeCand, check with the returned optimal ridge parameter. If the selected parameter is the maximum/minimum values in the candidates, it is possible that the selected one is too small/big.
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
#' sampWiener <- MakeFPCAInputs(IDs = rep(1:n, each=length(pts)), 
#'                              tVec = rep(pts, times = n), 
#'                              yVec = t(sampWiener))
#' res <- FOptDes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, p=2,
#'                isSequential=FALSE, RidgeCand = seq(0.02,0.2,0.02))
#' @references
#' \cite{Ji, H., Mueller, H.G. (2016) "Optimal Designs for Longitudinal and Functional Data" Journal of the Royal Statistical Society: Series B (Statistical Methodology)}
#' 
#' @export

FOptDes <- function(Ly, Lt, Resp, p = 3, optns = list(),
                        isRegression = !missing(Resp), isSequential = FALSE, RidgeCand = NULL){
  # check inputs
  if(is.null(RidgeCand)){
    stop("RidgeCand missing! Need to specify at least one ridge candidate.")
  }
  if( !(is.vector(RidgeCand) && is.numeric(RidgeCand)) ){
    stop("RidgeCand does not have the correct input format! Need to be a vector of positive numbers.")
  }
  if(any(RidgeCand <= 0)){
    stop("Some ridge candidates are non-positive! Change RidgeCand to make sure all ridge candidates are postive")
  }
  if( !(is.numeric(p) && p==as.integer(p) && p > 0) ){
    stop("Argument 'p' is not a positive integer! Need to specify a positive integer as the number of design points to be selected.")
  }
  if(isRegression){
    if(!(is.numeric(Resp) && is.vector(Resp)) ){
      stop("Resp does not have the correct input format! Need to be a vector of numbers.")
    }
    if( length(Resp) != length(Ly) ){
      stop("Resp does not have the same length as Ly! Double check the data inputs.")
    }
    cat("Finding optimal designs for scalar response prediction.\n")
  } else {
    cat("Finding optimal designs for trajectory recovery.\n")
  }
  
  CheckData(y = Ly, t = Lt);
  inputData <- HandleNumericsAndNAN(Ly, Lt);
  y <- inputData$Ly;
  t <- inputData$Lt;  
  
  obsGrid = sort(unique(c(unlist(t))));

  optns$nRegGrid = as.integer(1+diff(range(obsGrid))/min(diff(obsGrid))); 
  # make sure that FPCA workGrid is a denser grid of obsGrid for cv
  # if measurement times are random, bin data first before run the function.
  optns = SetOptions(y, t, optns);
  
  numOfCurves = length(y);
  CheckOptions(t, optns, numOfCurves);
  if(optns$dataType == "Dense"){
    isDense = TRUE;
  } else {
    isDense = FALSE; # currently dense with missing is treated as sparse
  }
  
  RegGrid = seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid);
  
  # find the best ridge parameter via cross validation
  if(length(RidgeCand) > 1){
    OptRidge <- MCVOptRidge(y = y, t = t, Resp = Resp, p = p, RidgeCand = RidgeCand,
                            isDense = isDense,
                            isRegression = isRegression, isSequential = isSequential)
    optridge <- OptRidge$optridge
  } else { # skip CV if only ridge is prespecified.
    cat("Only 1 ridge candidate in RidgeCand. The candidate is used and cross validation is skipped.\n")
    OptRidge <- RidgeCand
    optridge <- RidgeCand
  }
  # find optdes with optridge
  TrainFPCA <- FPCA(Ly=y, Lt=t, optns=optns)
  if(isRegression == FALSE){ # Trajectory Recovery
    BestDesTR <- BestDes_TR(p=p, ridge=optridge, workGrid=TrainFPCA$workGrid,
                            Cov=TrainFPCA$fittedCov, isSequential=isSequential)$best
    # calculate R2_X
    VarX <- sum(TrainFPCA$lambda)
    mu <- TrainFPCA$mu
    Cov <- TrainFPCA$fittedCov
    ridgeCov <- TrainFPCA$fittedCov + diag(optridge, nrow(Cov))
    R2XNum <- sum(diag(Cov[,BestDesTR] %*% solve(ridgeCov[BestDesTR, BestDesTR]) %*% Cov[BestDesTR,]))*diff(RegGrid)[1]
    R2X <- R2XNum/VarX
    R2Xadj <- 1-(1-R2X)*(length(y)-1)/(length(y)-p-1)
    if(R2X >= 1){
      warning("Coefficient of determination is greater than 1! Select other ridge candidates for proper regularization.")
    }
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
    BestDesSR <- BestDes_SR(p=p, ridge=optridge, workGrid=TrainFPCA$workGrid,
                            Cov=TrainFPCA$fittedCov, CCov=CCov, isSequential=isSequential)$best
    R2Y <- (var(Resp) - CCov[BestDesSR] %*% solve(ridgeCov[BestDesSR,BestDesSR]) %*% CCov[BestDesSR])/var(Resp)
    R2Yadj <- 1-(1-R2Y)*(length(y)-1)/(length(y)-p-1)
    if(R2Y >= 1){
      warning("Coefficient of determination is greater than 1! Select other ridge candidates for proper regularization.")
    }
    return(list(OptDes = RegGrid[BestDesSR], R2 = R2Y, R2adj = R2Yadj, OptRidge = OptRidge))
  }
}

