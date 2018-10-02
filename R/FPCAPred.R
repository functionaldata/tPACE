#'@title Functional Principal Components Score prediction

#'@description  This function performs prediction for new y and newt based on the returned fits from FPCA().

#'@param yfpca An object that is returned by FPCA().
#'@param Lnewy A length m list of new measurements for new subjects.
#'@param Lnewt A length m list of new time points for new subjects if all new subjects are evaluated at the same time, Lnewt can be a row vector of time points for one subject.
#'@param dataType 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' (see FPCA() for more details). Default is the same as in fpcaObj.

#'@return ypred: a length m list of predicted measurements for new subjects
#'xi_new: m x K matrix of new estimated FPC scores
#'xi_var: K*K matrix, Var(PC score)-Var(estimated PC score). The omega matrix in equation (7) of the paper, which is used to construct the point-wise C.I. for X_i[t].


FPCApred = function(fpcaObj,Lnewy,Lnewt,dataType){
  p = fpcaObj$optns
  if(is.na(regular)){
    dataType = optns$dataType
  }
  if(is.vector(Lnewt)){
    Lnewt = split(matrix(Lnewt,length(Lnewy),length(Lnewt),byrow = T),rep(1:length(Lnewt),each = length(Lnewt)))
  }
  mu = fpcaObj$mu
  phi = fpcaObj$phi
  lambda = fpcaObj$lambda
  sigma2 = fpcaObj$sigma2
  error = p$error
  method = p$methodXi
  shrink = p$shrink
  workgrid = fpcaObj$workGrid
  rho = fpcaObj$rho
  
  muSub = list()
  phiSub = list()
  k = ncol(phi)
  for(i in 1:length(Lnewt)){
    muSub = c(muSub,list(spline(x = workgrid,y=mu,xout = Lnewy[[i]])$y))
    phii = matrix(0,workgrid,k)
    for(j in 1:k){
      phii[,k] = spline(x = workgrid,y=phi[,k],xout = Lnewt[[i]])$y
    }
    phiSub = c(phiSub,list(phii))
  }
  
  ncohort = length(Lnewy)
  LAMBDA = diag(lambda)
  if(method == 'IN'){
    xi_var = c()
  }
  else{
    xi_var = list()
  }
  
  if(error == TRUE){
    y_predOrig = list()
    xi_est = matrix(0,ncohort,k)
    zeta_est = xi_est
    if(method == 'CE'){
      for(i in 1:ncohort){
        phii = phiSub[[i]]
        mui = muSub[[i]]
        yi = Lnewy[[i]]
        error0 = sigma2*diag(length(yi))
        A = (LAMBDA %*% t(phii)) %*% solve(phii %*% LAMBDA %*% t(phii)+error0)
        xi_est[i,] = (yi-mui) %*% A
        xi_var = c(xi_var,list(LAMBDA-A %*% (LAMBDA %*% t(phii))))
        y_predOrig = c(y_predOrig,list(mu_i+xi_est[i,] %*% phii))
      }
    }
    else{
        m = length(yi)
        for(i in 1:ncohort){
          phii = phiSub[[i]]
          mui = muSub[[i]]
          yi = Lnewy[[i]]
          for(j in 1:k){
            prod = (yi-mui) * phii[,j]
            if(shrink == 0){
              xi_est[i,j] = trapzRcpp(Lnewt[[i]],prod)
            }
            else{
              zeta_est[i,j] = trapzRcpp(Lnewt[[i]],prod)
              xi_est[i,j]=lambda[j]*zeta_est[i,j]/(lambda[k]+sigma2/m)
            }
          }
          y_predOrig = c(y_predOrig,list(mu_i+xi_est[i,] %*% phii))
        }
    }
  }
  else{
    y_predOrig = list()
    xi_est = matrix(0,ncohort,k)
    zeta_est = xi_est
    if(method == 'CE'){
      for(i in 1:ncohort){
        phii = phiSub[[i]]
        mui = muSub[[i]]
        yi = Lnewy[[i]]
        A = (LAMBDA %*% t(phii)) %*% solve(phii %*% LAMBDA %*% t(phii))
        xi_est[i,] = (yi-mui) %*% A
        xi_var = c(xi_var,list(LAMBDA-A %*% (LAMBDA %*% t(phii))))
        y_predOrig = c(y_predOrig,list(mu_i+xi_est[i,] %*% phii))
      }
    }
    else{
      for(i in 1:ncohort){
        for(j in 1:k){
          prod = (yi-mui) * phii[,j]
          xi_est[i,j] = trapzRcpp(Lnewt[[i]],prod)
        }
        y_predOrig = c(y_predOrig,list(mu_i+xi_est[i,] %*% phii))
      }
    }
  }
  return(list(ypred = y_predOrig,xi_est,xi_var))
}