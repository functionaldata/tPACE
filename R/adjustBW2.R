# This function adjusts the bandwidth for smoothing of covariance function in PCA
# bw_xcov is a vector input
##########################################################################
# Input:  - kernel: 'gauss', 'epan' or other kernel type
#         - bw_xcov: working bandwidth for covariance surface smoothing, a vector
#         - nder: degree of derivative need to be estimated
#         - regular: whether the functional data is dense and regular ("Dense")
#         - verbose: 
##########################################################################
# Output: - bw_xcov: bandwidth after adjustment, a vector
##########################################################################

adjustBW2 <- function(kernel, bw_xcov, nder, regular, verbose){
  # for Gaussian kernel
  if(kernel == 'gauss'){
    if(regular == "Dense"){
      bwxcov_fac = c(1.1, 0.8, 0.8)
    }
    else bwxcov_fac = c(1.1, 1.2, 2)
  
    if(nder > 2){ facID = 3 }
    else if(nder >= 0 && nder <= 2){facID = nder + 1}
    else facID = 1
  
  bw_xcov = bw_xcov * bwxcov_fac[facID]
  }
  # for Epanechnikov kernel
  if(kernel == 'epan'){
    if(regular == "Dense"){
      bwxcov_fac = c(1.1, 1.0, 1.1)
    }
    else bwxcov_fac = c(1.1, 1.2, 1.5)
    
    if(nder > 2){ facID = 3 }
    else if(nder >= 0 && nder <= 2){facID = nder + 1}
    else facID = 1
    
    bw_xcov = bw_xcov * bwxcov_fac[facID]    
  }
  bw_xcov
}  