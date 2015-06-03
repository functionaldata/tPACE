# This function adjusts the bandwidth for smoothing of covariance function in PCA
# bw_userCov is a vector input
##########################################################################
# Input:  - kernel: 'gauss', 'epan' or other kernel type
#         - bw_userCov: working bandwidth for covariance surface smoothing, a vector
#         - nder: degree of derivative need to be estimated
#         - dataType: whether the functional data is dense and dataType ("Dense")
#         - verbose: 
##########################################################################
# Output: - bw_userCov: bandwidth after adjustment, a vector
##########################################################################

adjustBW2 <- function(kernel, bw_userCov, nder, dataType, verbose){
  # for Gaussian kernel
  if(kernel == 'gauss'){
    if(dataType == "Dense"){
      bwuserCov_fac = c(1.1, 0.8, 0.8)
    }
    else bwuserCov_fac = c(1.1, 1.2, 2)
  
    if(nder > 2){ facID = 3 }
    else if(nder >= 0 && nder <= 2){facID = nder + 1}
    else facID = 1
  
  bw_userCov = bw_userCov * bwuserCov_fac[facID]
  }
  # for Epanechnikov kernel
  if(kernel == 'epan'){
    if(dataType == "Dense"){
      bwuserCov_fac = c(1.1, 1.0, 1.1)
    }
    else bwuserCov_fac = c(1.1, 1.2, 1.5)
    
    if(nder > 2){ facID = 3 }
    else if(nder >= 0 && nder <= 2){facID = nder + 1}
    else facID = 1
    
    bw_userCov = bw_userCov * bwuserCov_fac[facID]    
  }
  bw_userCov
}  