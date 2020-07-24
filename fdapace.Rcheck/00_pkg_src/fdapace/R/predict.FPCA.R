#' Predict FPC scores and curves for a new sample given an FPCA object
#'
#' Return a list containing the matrix with the first k FPC scores according to conditional expectation or numerical integration, the matrix of predicted trajectories and the prediction work grid.
#'
#' @param object An FPCA object.
#' @param newLy  A list of \emph{n} vectors containing the observed values for each individual.
#' @param newLt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param sigma2 The user-defined measurement error variance. A positive scalar. (default: rho if applicable, otherwise sigma2 if applicable, otherwise 0 if no error. )
#' @param K The scalar defining the number of clusters to define; (default: 1).
#' @param xiMethod The integration method used to calculate the functional principal component scores 
#' (standard numerical integration 'IN' or conditional expectation 'CE'); default: 'CE'.
#' @param ... Not used.
#' 
#' @return  A list containing the following fields:
#' \item{scores}{A matrix of size \emph{n}-by-\emph{k} which comprise of the predicted functional principal component scores.}
#' \item{predCurves}{A matrix of size \emph{n}-by-\emph{l} where \emph{l} is the length of the work grid in \emph{object}.}
#' \item{predGrid}{A vector of length \emph{l} which is the output grid of the predicted curves. This is same is the workgrid of \emph{object}.}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50
#' pts <- seq(0, 1, by=0.05)
#' # The first n samples are for training and the rest testing
#' sampWiener <- Wiener(2 * n, pts)
#' sparsity <- 2:5
#' train <- Sparsify(sampWiener[seq_len(n), , drop=FALSE], pts, sparsity)
#' test <- Sparsify(sampWiener[seq(n + 1, 2 * n), , drop=FALSE], pts, sparsity)
#' res <- FPCA(train$Ly, train$Lt)
#' pred <- predict(res, test$Ly, test$Lt, K=5)
#' plot(pred$predGrid, pred$predCurves[1,])
#' }
#' 
#' @method predict FPCA
#' @export

predict.FPCA <- function(object, newLy, newLt, sigma2 = NULL, K = 1, xiMethod = 'CE', ...){
  fpcaObj = object;
  
  if(class(fpcaObj) != "FPCA"){
    stop('Please provide a valid FPCA object.')
  } 
  
  if(! all.equal( sapply(newLt, length), sapply(newLy, length) ) ){
    stop('The size of the vectors in newLt and newLy differ. They must be equal.')
  } 
  
  # Standard data checks/massasing as in FPCA()
  CheckData(newLy, newLt)
  inputData <- HandleNumericsAndNAN(Lt = newLt, Ly = newLy);
  newLy <- inputData$Ly;
  newLt <- inputData$Lt;
  
  if(is.null(sigma2)){
    sigma2 <- ifelse( !is.null(fpcaObj$rho), fpcaObj$rho, 
                      ifelse( !is.null(fpcaObj$sigma2), fpcaObj$sigma2,
                              ifelse( !fpcaObj$optns$error, 0, stop('sigma2 cannot be determined.'))))  
  } else {
    if(!is.numeric(sigma2) ){
      stop('sigma2 is not numeric, we sure about this? :D')
    }
  }
  
  if(K > fpcaObj$selectK){
    stop( paste0( collapse = '', 'You cannot get FPC scores for more components than what is already available. (', fpcaObj$selectK ,').' ))
  } 
  
  if( !(xiMethod %in% c('CE','IN')) ){
    stop( paste0( collapse = '', 'Unrecognised method to calculate the FPC scores.'))
  }
  if( (fpcaObj$optns$dataType == 'Sparse') && (xiMethod == 'IN')){
    stop( 'Trapezoid Numerical intergration (IN) is invalid for sparse data.')
  }
  
  rangeOfOrigData <- range(unlist(fpcaObj$workGrid))
  rangeOfNewData <- range(unlist(newLt))
  
  if(rangeOfNewData[1] < rangeOfOrigData[1]){
    stop("The new data's lower range is below the original data's lower range.")
  } else{
    if(rangeOfNewData[2] > rangeOfOrigData[2]){
      stop("The new data's upper range is above the original data's upper range.")
    }
  }
  
  MuObs = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = fpcaObj$obsGrid, mu = fpcaObj$mu)
  CovObs = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = fpcaObj$obsGrid, Cov =  fpcaObj$fittedCov)
  PhiObs = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = fpcaObj$obsGrid, phi = fpcaObj$phi)
  
  # Get scores  
  if ( xiMethod == 'CE') {
    scoresObj <- GetCEScores(y = newLy, t = newLt, optns = fpcaObj$optns, mu = MuObs, obsGrid = fpcaObj$obsGrid, sigma2 = sigma2,
                             fittedCov = CovObs, lambda = fpcaObj$lambda, phi = PhiObs)
    finalXiEst <- t(do.call(cbind, scoresObj['xiEst', ]))[, seq_len(K), drop=FALSE]
  } else if (xiMethod == 'IN') {
    scoresObj <- mapply(function(yvec,tvec)
      GetINScores(yvec, tvec,optns = fpcaObj$optns,obsGrid = fpcaObj$obsGrid,mu = MuObs,lambda =fpcaObj$lambda ,phi = PhiObs,sigma2 = sigma2),newLy,newLt)
    finalXiEst <- t(do.call(cbind, scoresObj['xiEst', ]))[, seq_len(K), drop=FALSE]
  }
  
  #return(finalXiEst)
  
  #Get predicted trajectories
  trajectoryEst <- t(fpcaObj$mu+fpcaObj$phi[,1:K]%*%t(finalXiEst))
  
  return(list(scores = finalXiEst, predCurves = trajectoryEst, predGrid = fpcaObj$workGrid))
  
}
