#' Create the options list used by FPCA
#' 
#' @param bwcov : bandwidth value for covariance function; positive numeric - default: determine automatically based on 'bwcovMethod'
#' @param bwcovMethod : bandwidth choice method for covariance function; 'GMeanAndGCV','CV','GCV - default: 'GMeanAndGCV'')
#' @param bwmu : bandwidth value for mean function is using CV or GCV; positive numeric - default: determine automatically based on 'bwmuMethod'
#' @param bwmuMethod : bandwidth choice method for mean function; 'GMeanAndGCV','CV','GCV - default: 'GMeanAndGCV''
#' @param dataType : do we have sparse or dense functional data; 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' - default:  determine automatically based on 'IsRegular' 
# '@param diagnosticsPlot : make diagnostics plot (design plot, mean, scree plot and first k (<=3) eigenfunctions); logical - default: FALSE}
#' @param error : assume measurement error in the dataset; logical - default: TRUE
#' @param FVEthreshold : Fraction-of-Variance-Explained threshold used during the SVD of the fitted covar. function; numeric (0,1] - default: 0.9999
#' @param kernel : smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "epan" for dense data else "gauss"
#' @param methodCov :  method to estimate covariance; 'PACE','RARE','CrossSectional' - automatically determined, user input ignored
#' @param methodMu :  method to estimate mu; 'PACE','RARE','CrossSectional' - automatically determined, user input ignored 
#' @param maxK : maximum number of principal components to consider; positive integer - default: min(20, N-1), N : # of curves
#' @param methodXi : method to estimate the PC scores; 'CE', 'IN' - default: 'CE'
#' @param ntest1 : number of curves used for CV when choosing bandwidth; [1,N] - default: min(30, N-1), N : # of curves
#' @param nRegGrid : number of support points in each direction of covariance surface; numeric - default: 51
#' @param numBins : number of bins to bin the data into; positive integer > 10, default: NULL 
#' @param selectionMethod : the method of choosing the number of principal components K; 'FVE','AIC','BIC' : default 'FVE' - only 'FVE' avaiable now/ default 'FVE')
#' @param shrink : apply shrinkage to estimates of random coefficients (dense data only); logical - default: FALSE
#' @param outPercent : 2-element vector in [0,1] indicating the outPercent data in the boundary - default (0,1)
#' @param rho : truncation threshold for the iterative residual. 'cv': choose rho by leave-one-observation out cross-validation; 'no': use the iterative sigma2 estimate - default "cv".
#' @param rotationCut : 2-element vector in [0,1] indicating the percent of data truncated during sigma^2 estimation; default  (0.25, 0.75))
#' @param useBinnedData : 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin) - default: 'AUTO'
#' @param useBins: testing purpose: whether to bin the same observed time points when 2D smoothing; logical - default: FALSE
#' @param userCov : user-defined smoothed covariance function; numerical matrix - default: NULL
#' @param userMu : user-defined smoothed mean function; numerical vector - default: NULL
#' @param verbose : display diagnostic messages; logical - default: FALSE
#' @return an option list
#' @examples 
#' optLst = CreateOptions(kernel='rect');  # Create options list with rectangular kernel 



CreateOptions = function(bwmu = NULL, bwmuMethod = NULL, bwuserCov = NULL, bwuserCovGcv = NULL, 
    ntest1 = NULL,  selectionMethod = NULL, FVEthreshold = NULL, numComponents = NULL,
    maxK = NULL, dataType = NULL, error = NULL, nRegGrid = NULL,
    methodXi = NULL, shrink = NULL,# newdata = NULL,
    kernel = NULL, 
    numBins = NULL, yname = NULL, # screePlot = NULL, designPlot = NULL, corrPlot = NULL, corrPlotType =NULL,
    diagnosticsPlot = NULL,
      rho = NULL, verbose = NULL, userMu = NULL, userCov = NULL, methodCov = NULL,
    methodMu = NULL, outPercent = NULL, useBinnedData = NULL, rotationCut = NULL){ 

 return( list(bwmu = bwmu, bwmuMethod = bwmuMethod, bwuserCov = bwuserCov, bwuserCovGcv = bwuserCovGcv,
          ntest1 = ntest1,  selectionMethod = selectionMethod, FVEthreshold = FVEthreshold,
          maxK = maxK, dataType = dataType, error = error, nRegGrid = nRegGrid, 
          methodXi = methodXi, shrink = shrink,# newdata = newdata,
          kernel = kernel, 
          # corrPlot = corrPlot, corrPlotType = corrPlotType, screePlot = screePlot, designPlot = designPlot, 
          numComponents = numComponents, diagnosticsPlot = diagnosticsPlot,
          numBins = numBins, yname = yname, rho = rho, rotationCut = rotationCut,
          verbose = verbose, userMu = userMu, userCov = userCov, methodMu= methodMu,  methodCov= methodCov, 
          outPercent = outPercent, useBinnedData = useBinnedData, rotationCut = rotationCut) )
}
