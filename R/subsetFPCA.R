### subset the FPCA object with specified number of components k
### and return the subsetted fpcaObj

subsetFPCA <- function(fpcaObj, k){
  fpcaObj$lambda <- fpcaObj$lambda[1:k]
  fpcaObj$phi <- fpcaObj$phi[,1:k, drop=FALSE]
  fpcaObj$xiEst <- fpcaObj$xiEst[,1:k, drop=FALSE]
  return(fpcaObj)
}