#' Selects number of functional principal components for
#' given FPCA output and selection criteria
#'
#' @param fpcaObj A list containing FPCA related subjects returned by MakeFPCAResults().
#' @param criterion A string or positive integer specifying selection criterion for number of functional principal components, available options: 'FVE', 'AIC', 'BIC', or the specified number of components - default: 'FVE'
#' @param FVEthreshold A threshold percentage specified by user when using "FVE" as selection criterion: (0,1] - default: NULL
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual - default: NULL
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to Ly - default: NULL
#'
#' @return A list including the following two fields:
#' \item{K}{An integer indicating the selected number of components based on given criterion.}
#' \item{criterion}{The calculated criterion value for the selected number of components, i.e. FVE, AIC or BIC value, NULL for fixedK criterion.}
#' \item{k}{Same as \code{K} for compatibility. WARNING: This will be removed in the next iteration}
#'
#' @export

SelectK = function(fpcaObj, criterion = 'FVE', FVEthreshold = 0.95, Ly = NULL, Lt = NULL){
  if(class(fpcaObj) != 'FPCA'){
    stop('Invalid Input: not a FPCA object!')
  }
  if(is.null(criterion)){
    stop('Invalid selection criterion. Selection criterion must not be NULL!')
  }
  if (length(criterion) != 1) {
    stop('The length of criterion needs to be 1')
  }
  if(!(criterion %in% c('FVE', 'AIC', 'BIC'))){
    if(is.numeric(criterion)){
      if(as.integer(criterion) != criterion || criterion <= 0){
        stop('Invalid selection criterion. To select fixed number of component, criterion needs to be a positive integer.')
      }
    } else {
      stop('Invalid selection criterion. Need to be one of "FVE", "AIC", "BIC" or a positive integer!')
    }
  }
  
  if(criterion %in% c('AIC','BIC')) {
    if(fpcaObj$optns$lean == TRUE && (is.null(Ly) || is.null(Lt))){
    stop("Option lean is TRUE, need input data Ly and measurement time list Lt to calculate log-likelihood.")
    }
    if(fpcaObj$optns$lean == FALSE){
      Ly <- fpcaObj$inputData$Ly
      Lt <- fpcaObj$inputData$Lt
    }
    if(criterion == 'AIC'){C = 2}
    else {C = log(length(Ly))}
    # FVE is not the selection criterion
    IC = rep(Inf, length(fpcaObj$lambda))
    for(i in 1:length(fpcaObj$lambda)){
      logliktemp = GetLogLik(fpcaObj, i, Ly = Ly, Lt = Lt)
      if(is.null(logliktemp)){
        stop('The covariance matrix of the estimated function is nearly singular! AIC or BIC is not applicable.')
      }
      IC[i] = logliktemp + C * i
      if(i > 1 && IC[i] > IC[i-1]){
        # cease whenever AIC/BIC stops decreasing
        K <- i-1
        criterion <- IC[i-1]
      } else if(i == length(fpcaObj$lambda)){
        K <- i
        criterion <- IC[i]
      }
    }
    #if(criterion != 'FVE'){
    #  return(K = length(fpcaObj$lambda))
    #}
  } else if(criterion == 'FVE'){
    # select FVE based on cumFVE in fpcaObj and specified FVEthreshold
    if(is.null(FVEthreshold)){stop('Need to specify FVEthreshold to choose number of components via FVE.')}
    cumFVE = fpcaObj$cumFVE
    buff <- .Machine[['double.eps']] * 100
    K <- min( which(cumFVE > FVEthreshold * 100 - buff) )
    criterion <- cumFVE[min(which(cumFVE > FVEthreshold * 100 - buff))]
  } else if (is.numeric(criterion) && criterion > 0) { # fixed K is specified.
    if(criterion > length(fpcaObj$lambda)){
      stop("Specified number of components is more than available components.")
    }
    K <- criterion
    criterion <- NULL
  } else {
    stop('Unknown criterion!')
  }
  
  # For compatibility reason, k is also returned.
  return(list(K=K, criterion=criterion, k=K)) 
}
