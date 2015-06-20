# This function converts dense regular functional input list
# to a matrix for easier dense case implementation
##########################################################################
# Input:  - y: list of n dense regular observed p-dim functional objects
##########################################################################
# Output: - ymat: n by p matrix containing all functional data
##########################################################################

List2Mat <- function(y){
  n = length(y)
  ymat = matrix(unlist(y), nrow = n, byrow = TRUE)
  
  return(ymat)
}