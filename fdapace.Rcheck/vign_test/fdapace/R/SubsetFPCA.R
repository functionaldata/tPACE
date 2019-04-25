### subset the FPCA object with specified number of components K
### and return the subsetted fpcaObj

SubsetFPCA <- function(fpcaObj, K){
  fpcaObj$lambda <- fpcaObj$lambda[1:K]
  fpcaObj$phi <- fpcaObj$phi[,1:K, drop=FALSE]
  fpcaObj$xiEst <- fpcaObj$xiEst[,1:K, drop=FALSE]
  fpcaObj$FVE <- fpcaObj$cumFVE[K]
  fpcaObj$cumFVE <- fpcaObj$cumFVE[1:K]
  return(fpcaObj)
}
