#' Predict FPC scores for a new sample given an FPCA object
#'
#' Return a matrix with the first k FPC scores according to conditional expectation or numerical integration.
#'
#' @param object An FPCA object.
#' @param newY  A list of \emph{n} vectors containing the observed values for each individual.
#' @param newX A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param sigma2 The user-defined measurement error variance. A positive scalar. (default: rho if applicable, otherwise sigma2 if applicable, otherwise 0 if no error. )
#' @param k The scalar defining the number of clusters to define; (default: 1).
#' @param im The integration method used to calculate the functional principal component scores ( standard numerical integration 'IN' or conditional expectation 'CE'); default: 'CE'.
#' @param ... Not used.
#' 
#' @return  A matrix of size n-by-k 
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt)
#' res
#'
#' @method predict FPCA
#' @export

predict.FPCA <- function(object, newX, newY, sigma2 = NULL, k = 1, im = 'CE', ...){
  fpcaObj = object;
  
  if(class(fpcaObj) != "FPCA"){
    stop('Please provide a valid FPCA object.')
  } 
  
  if(! all.equal( sapply(newX, length), sapply(newY, length) ) ){
    stop('The size of the vectors in newX and newY differ. They must be equal.')
  } 
  
  # Standard data checks/massasing as in FPCA()
  CheckData(newY, newX)
  inputData <- HandleNumericsAndNAN(Lt = newX, Ly = newY);
  newY <- inputData$Ly;
  newX <- inputData$Lt;
  
  if(is.null(sigma2)){
    sigma2 <- ifelse( !is.null(fpcaObj$rho), fpcaObj$rho, 
                      ifelse( !is.null(fpcaObj$sigma2), fpcaObj$sigma2,
                              ifelse( !fpcaObj$optns$error, 0, stop('sigma2 cannot be determined.'))))  
  } else {
    if(!is.numeric(sigma2) ){
      stop('sigma2 is not numeric, we sure about this? :D')
    }
  }
  
  if(k > fpcaObj$selectK){
    stop( paste0( collapse = '', 'You cannot get FPC scores for more components than what is already available. (', fpcaObj$selectK ,').' ))
  } 
  
  if( !(im %in% c('CE','IN')) ){
    stop( paste0( collapse = '', 'Unrecognised method to calculate the FPC scores.'))
  }
  if( (fpcaObj$optns$dataType == 'Sparse') && (im == 'IN')){
    stop( 'Trapezoid Numerical intergration (IN) is invalid for sparse data.')
  }
  
  rangeOfOrigData <- range(unlist(fpcaObj$workGrid))
  rangeOfNewData <- range(unlist(newX))
  
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
  if ( im == 'CE') {
    scoresObj <- GetCEScores(y = newY, t = newX, optns = fpcaObj$optns, mu = MuObs, obsGrid = fpcaObj$obsGrid, sigma2 = sigma2,
                             fittedCov = CovObs, lambda = fpcaObj$lambda, phi = PhiObs)
  finalXiEst <- t(do.call(cbind, scoresObj['xiEst', ]))[,1:k] 
  } else if (im == 'IN') {
    ymat = List2Mat(newY,newX)
    scoresObj <- GetINScores(ymat = ymat, t = newX,optns = fpcaObj$optns,mu = MuObs,lambda =fpcaObj$lambda ,phi = PhiObs,sigma2 = sigma2)
    finalXiEst <- scoresObj$xiEst[,1:k]
  }
  
  return(finalXiEst)
}
