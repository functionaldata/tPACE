#' Create the PCA option list
#' 
#' @param bwmu : bandwidth choice for mean function is using CV or GCV
#' @param bwmu_gcv : bandwidth choice method for mean function is GCV if bwmu = 0
#' @param bwcov : bandwidth choice for covariance function is CV or GCV
#' @param bwcov_gcv : bandwidth choice method for covariance function is GCV if bwxcov = c(NULL,0)
#' @param ntest1 : number of curves used for CV when choosing bandwidth 
#' @param ngrid1 : number of support points for the covariance surface 
#' @param selection_k : the method of choosing the number of principal components K
#' @param FVE_threshold : Fraction-of-Variance-Explained
#' @param maxk : maximum number of principal components to consider
#' @param regular : do we have regular or sparse functional data
#' @param error : error assumption with measurement error
#' @param ngrid : number of support points in each direction of covariance surface 
#' @param method : method to estimate the PC scores
#' @param shrink : apply shrinkage to estimates of random coefficients (regular data only)
#' @param newdata : new data points to estimate
#' @param kernel : smoothing kernel choice
#' @param numBins : number of bins
#' @param yname : name of the variable analysed
#' @param screePlot : make scree plot
#' @param designPlot : make design plot
#' @param corrPlot : make correlation plot
#' @param rho : truncation threshold for the iterative residual  
#' @param verbose : display diagnostic messages
#' @param xmu : user-defined smoothed mean function 
#' @param xcov : user-defined smoothed covariance function
#' @param method_mu :  method to estimate mu
#' @param out_percent : number in [NULL,1] indicating the out_percent data in the boundary
#' @param use_binned_data : 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin)
#' @return an option list
#' @examples 
#' 1 + 3

CreateOptions = function(bwmu = NULL, bwmu_gcv = NULL, bwxcov = NULL, bwxcov_gcv = NULL, 
    ntest1 = NULL, ngrid1 = NULL, selection_k = NULL, FVE_threshold = NULL,
    maxk = NULL, regular = NULL, error = NULL, ngrid = NULL,
    method = NULL, shrink = NULL, newdata = NULL, kernel = NULL, 
    numBins = NULL, yname = NULL, screePlot = NULL, designPlot = NULL, 
    corrPlot = NULL,     rho = NULL, verbose = NULL, xmu = NULL, xcov = NULL, 
    method_mu = NULL, out_percent = NULL, use_binned_data = NULL){ 

 return( list(bwmu = bwmu, bwmu_gcv = bwmu_gcv, bwxcov = bwxcov, bwxcov_gcv = bwxcov_gcv,
          ntest1 = ntest1, ngrid1 = ngrid1, selection_k = selection_k, FVE_threshold = FVE_threshold,
          maxk = maxk, regular = regular, error = error, ngrid = ngrid, 
          method = method, shrink = shrink, newdata = newdata, kernel = kernel, corrPlot = corrPlot,	
          numBins = numBins, yname = yname, screePlot = screePlot, designPlot = designPlot, rho = rho, 
          verbose = verbose, xmu = xmu, xcov = xcov, method_mu= method_mu, out_percent = out_percent, use_binned_data = use_binned_data) )
}
