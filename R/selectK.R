#' Selects number of functional principal components for
#' given FPCA output and selection criteria
#'
#' @param ret A list containing FPCA related subjects returned by MakeFPCAResults().
#' @param criterion A string specifying selection criterion for number of functional principal components, available options: 'FVE', 'AIC', 'BIC', 'fixedK' - default: 'AIC'
#' @param FVEthreshold A threshold percentage specified by user when using "FVE" as selection criterion - default: NULL
#' @param fixedK An integer: user-specified number of components to be chosen - default: NULL
#' @param y A list of \emph{n} vectors containing the observed values for each individual - default: NULL
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y - default: NULL
#'
#' @return A list including the following two fields:
#' \item{k}{An integer indicating the selected number of components based on given criterion.}
#' \item{criterion}{The calculated criterion value for the selected number of components, i.e. FVE, AIC or BIC value, NULL for fixedK criterion.}
#'
#' @export

selectK = function(ret, criterion = 'AIC', FVEthreshold = NULL, fixedK = NULL, y = NULL, t = NULL){
  if(class(ret) != 'FPCA'){stop('Invalid Input: not a FPCA object!')}
  if(!(criterion %in% c('FVE', 'AIC', 'BIC', 'fixedK') || is.null(criterion))){
    stop('Invalid selection criterion. Need to be one of "FVE", "AIC", "BIC" or "fixedK"!')
  }
  # if called within FPCA function
  if(is.null(criterion)){
    if(ret$optns$selectionMethod %in% c('AIC','BIC')) {
      if(ret$optns$lean == TRUE && (is.null(y) || is.null(t))){
        stop("Option lean is TRUE, need input data y and measurement time list t to calculate log-likelihood.")
      }
      if(ret$optns$lean == FALSE){
        y <- ret$inputData$y
        t <- ret$inputData$t
      }
      if(ret$optns$selectionMethod == 'AIC'){C = 2}
      else {C = log(length(y))}
      # FVE is not the selection criterion
      IC = rep(Inf, length(ret$lambda))
      for(i in 1:length(ret$lambda)){
        logliktemp = getLogLik(ret, i, y = y, t = t)
        if(is.null(logliktemp)){
          stop('The covariance matrix of the estimated function is nearly singular! AIC or BIC is not applicable.')
        }
        IC[i] = logliktemp + C * i
        if(i > 1 && IC[i] > IC[i-1]){
          # cease whenever AIC/BIC stops decreasing
          return(list(k = i-1, criterion = IC[i-1]))
        }
        if(i == length(ret$lambda)){
          return(list(k = i, criterion = IC[i]))
        }
      }
      #if(criterion != 'FVE'){
      #  return(k = length(ret$lambda))
      #}
    } else if(ret$optns$selectionMethod == 'FVE'){
      # no selection based on FVE, already selected in GetEigenAnalysisResults
      return( list(k = length(ret$lambda), criterion = ret$cumFVE[length(ret$lambda)]))
    } else { # fixed K is specified.
      if(is.null(ret$optns$fixedK)){stop("Need to specify fixedK in optns when selectionMethod is 'fixedK'.")}
      if(fixedK > length(ret$lambda)){
        stop("Specified number of components fixedK is more than available components.")
      }
      return(list(k = fixedK, criterion = NULL))
    }
  } else {  # if called outside FPCA by user
  if(criterion %in% c('AIC','BIC')) {
    if(ret$optns$lean == TRUE && (is.null(y) || is.null(t))){
      stop("Option lean is TRUE, need input data y and measurement time list t to calculate log-likelihood.")
    }
    if(ret$optns$lean == FALSE){
      y <- ret$inputData$y
      t <- ret$inputData$t
    }
    if(criterion == 'AIC'){C = 2}
    else {C = log(length(y))}
    # FVE is not the selection criterion
    IC = rep(Inf, length(ret$lambda))
    for(i in 1:length(ret$lambda)){
      logliktemp = getLogLik(ret, i, y = y, t = t)
      if(is.null(logliktemp)){
        stop('The covariance matrix of the estimated function is nearly singular! AIC or BIC is not applicable.')
      }
      IC[i] = logliktemp + C * i
      if(i > 1 && IC[i] > IC[i-1]){
        # cease whenever AIC/BIC stops decreasing
        return(list(k = i-1, criterion = IC[i-1]))
      }
      if(i == length(ret$lambda)){
        return(list(k = i, criterion = IC[i]))
      }
    }
    #if(criterion != 'FVE'){
    #  return(k = length(ret$lambda))
    #}
  } else if(criterion == 'FVE'){
    # FVE is the selection criterion
    if(is.null(FVEthreshold)){stop('Need to specify FVE threshold to choose number of components.')}
    cumFVE = ret$cumFVE
    return( list(k = min( which(cumFVE > FVEthreshold * 100) ), criterion = cumFVE[min(which(cumFVE > FVEthreshold))]))
  } else { # fixed K is specified.
    if(is.null(fixedK)){stop("Need to provide fixedK as function inpit when criterion is 'fixedK'.")}
    if(fixedK > length(ret$lambda)){
      stop("Specified number of components fixedK is more than available components.")
      fixedK = length(ret$lambda)
    }
    return(list(k = fixedK, criterion = NULL))
  }
  }
}
