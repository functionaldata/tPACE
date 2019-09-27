#'@title Conditional Quantile estimation with functional covariates
#'
#'@description  Main function to implement conditional Quantile estimation with functional covariates and scalar response. The method includes 3 steps: 
#'1) FPCA using the PACE method for X(t_x)
#'2) Computation of the conditional distribution function through a functional generalized linear model.
#'3) Prediction of quantiles for given predictor values
#'
#'@param Lx A length n list of predictor function where x[[i]] is the row vector of measurements for ith subject, i=1,...,n
#'@param Lt_x A length n list where the observations of x are taken, t_x[[i]] is a row vector of time points where x[[i]] are observed, i=1,...,n
#'@param y A 1*n vector for scalar response y. y[i] is the response value for the ith subject, i = 1,...,n.
#'@param outQ A vector of desired quantile levels with default value outQ = c(0.1, 0.25, 0.5, 0.75, 0.9).
#'@param optns_x A list of options for predictor x with control parameters specified by list(name=value) with default NA. See function FPCA for details.
#'@param isNewSub A 1*n vector of 0s or 1s, where n is the total count of subjects. 0 denotes the corresponding subject is only used for training and 1 denotes the corresponding subject is only used for prediction. (default: 0's)
#'
#'@return A list of the following
#'\item{pred_quantile}{A matrix of n*length(outQ) where the the first nn (number of 0s in \code{isNewSub}) rows containing fitted  conditional quantiles of Y corresponding to the training subjects, and the last n-nn rows containing predicted conditional quantiles of Y corresponding to the subjects isNewSub ==1.}
#'\item{pred_CDF}{A matrix of n*100. The ith row contains the fitted or predicted conditional distribution function \eqn{F(y|X_i)}, evaluated at an equally spaced grid of 100 points.}
#'\item{b}{A matrix of 50*(K+1) contains the coefficient functions, defined as \eqn{F(y|X) = g(\sum_(k=0)^K b_k(y)\xi_k)}, see equation (5) in the paper for details, where K is the number of components selected to expand the predictor functions X, and \eqn{\xi_k} is the kth principal component score.}
#'
#'@examples
#'set.seed(10)
#'
#'n = 200
#'npred = 50
#'m = 50
#'xi <- Wiener(n, 0:m/m)
#'
#'x=list()
#'t_x=list()
#'y=numeric(n)
#'for(i in 1:n){
#'  t_x = c(t_x,list(0:m/m))
#'  x = c(x,list(xi[i,]))
#'  y[i] = 5*rnorm(1)+2*sum(xi[i,])
#'}
#'
#'outQ = c(0.1,0.25,0.5,0.75,0.9,0.95)
#'isNewSub = c(rep(0,150),rep(1,50))
#'qtreg = FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewSub)
#'
#'@references
#' \cite{Chen, K., Müller, H.G. (2011). Conditional quantile analysis when covariates are functions, with application to growth data. 
#' J. Royal Statistical Society B 74, 67-89}
#' @export


FPCquantile = function(Lx,Lt_x,y,outQ=c(0.1,0.25,0.5,0.75,0.9),optns_x = NULL,isNewSub = NULL){
  CheckData(Lx,Lt_x)
  if(min(outQ) < 0 || max(outQ) > 1){
    stop("quantile levels must be between 0 and 1")
  }
  if(length(Lx) != length(y)){
    stop('length of Lx and y must agree')
  }
  n = length(Lx)
  if(is.null(isNewSub)){
    isNewSub = rep(0,n)
  }
  nn = sum(isNewSub == 0)
  Lx_e = list()
  Lt_x_e = list()
  Lnewx = list()
  Lnewt_x = list()
  for(i in 1:n){
    if(isNewSub[i] == 0){
      Lx_e = c(Lx_e,list(Lx[[i]]))
      Lt_x_e = c(Lt_x_e,list(Lt_x[[i]]))
    }else{
      Lnewx = c(Lnewx,list(Lx[[i]]))
      Lnewt_x = c(Lnewt_x,list(Lt_x[[i]]))
    }
  }
  ty= y[isNewSub==0]
  optns_x = SetOptions(Lx, Lt_x, optns_x)
  CheckOptions(Lt_x, optns_x,n)
  verbose_x = optns_x$verbose
  if(verbose_x == TRUE){
    print('Obtain functional object for x:\n')
  }
  xfpca = FPCA(Lx_e, Lt_x_e, optns_x) 
  xitrain = xfpca$xiEst
  optns = xfpca$optns

  zm = 50
  ysort = sort(ty)
  zgrid= seq(ysort[5], ysort[nn-5], length.out=zm)
  K_x = ncol(xitrain)
  b = matrix(0,zm,K_x+1)
  for(i in 1:zm){
    z=zgrid[i]
    ytrain =array(0,nn)
    ytrain[ty<=z]=1
    glmfit= glm(ytrain ~ xitrain, family = binomial(link= 'logit'))
    b[i,] = glmfit$coefficients
  }
  xi = matrix(0,n, K_x)
  xi[1:nn,] = xitrain
  if(length(Lnewx) > 0){
    xi[(nn+1):n,] = predict.FPCA(xfpca,Lnewx,Lnewt_x,K = K_x,xiMethod = optns$methodXi)
  }
  
  m = 100
  outF = seq(ysort[5], ysort[nn-5], length.out=m)
  predF = matrix(0,n, m)
  nQ = length(outQ)
  predQ = matrix(0,n, nQ)
  
  bw = 4*(zgrid[2]-zgrid[1])
  for(i in 1:n){
    CDF = 1:zm
    for(j in 1:zm){
      eta = b[j,] %*% c(1,xi[i,])
      if(eta>=30){
        CDF[j] = 1
      }else if(eta <= -30){
        CDF[j] = 0
      }else{
        CDF[j] = exp(eta)/(1+exp(eta))
      }
    }
    win=array(1,zm);
    predF[i,] =Lwls1D(bw,'epan', win, zgrid, CDF, outF)
  }
  
  for(i in 1:n){
    for(j in 1:nQ){
      if(predF[i,m] <= outQ[j]){
        predQ[i,j]= outF[m]
      }else{
        predQ[i,j]= outF[which(predF[i,] > outQ[j])[1]]
      }
    }
  }
  return(list(pred_quantile = predQ,pred_CDF = predF,beta = b))
}
