#'@title Conditional Quantile estimation with functional covariates

#'@description  Main function to implement conditional Quantile estimation with functional covariates and scalar response. The method includes 3 steps: 
#'1) FPCA using the PACE method for X(t_x)
#'2) Computation of the conditional distribution function through a functional generalized linear model.
#'3) Prediction of quantiles for given predictor values

#'@param Lx A length n list of predictor function where x[[i]] is the row vector of measurements for ith subject, i=1,...,n
#'@param Lt_x A length n list where the observations of x are taken, t_x[[i]] is a row vector of time points where x[[i]] are observed, i=1,...,n
#'@param y A 1*n vector for scalar response y. y[i] is the response value for the ith subject, i = 1,...,n.
#'@param outQ A vector of desired quantile levels with default value outQ = c(0.1, 0.25, 0.5, 0.75, 0.9).
#'@param optns_x A list of options for predictor x with control parameters specified by list(name=value) with default NA. See function FPCA for details.
#'@param isNewSub A 1*n vector of 0s or 1s, where n is the total count of subjects. 0 denotes the corresponding subject is only used for estimation and 1 denotes the corresponding subject is only used for prediction. (default: 0's)

#'@return predQ: a matrix of n*length(outQ) where the the first nn (number of 0s in isNewSub) rows containing fitted  conditional quantiles of Y corresponding to the trainning subjects, and the last n-nn rows containing predicted conditional quantiles of Y corresponding to the subjects isNewSub ==1.
#'predF: a matrix of n*100. The ith row contains the fitted or predicted conditional distribution function F(y|X_i), evaluated at an equally spaced grid of 100 points.
#'b: a matrix of 50*(K+1) contains the coefficient functions, defined as F(y|X) = g(\sum_(k=0)^K b_k(y)\xi_k), see equation (5) in the paper for details, where K is the number of components selected to expand the predictor functions X, and \xi_k is the kth principal component score.

#'@references Chen, K., M\"uller, H.G. (2011). Conditional quantile analysis when covariates are functions, with application to growth data. J. Royal Statistical Society B.


FPCquantile = function(Lx,Lt_x,y,outQ=c(0.1,0.25,0.5,0.75,0.9),optns_x,isNewsub){
  CheckData(Lx,Lt_x)
  n = length(Lx)
  if(is.na(isNewSub)){
    isNewsub = rep(0,n)
  }
  nn = sum(isNewsub == 0)
  Lx_e = list()
  Lt_x_e = list()
  Lnewx = list()
  Lnewt_x = list()
  for(i in 1:n){
    if(isNewsub == 0){
      Lx_e = c(Lx_e,list(Lx[[i]]))
      Lt_x_e = c(Lt_x_e,list(Lt_x[[i]]))
    }
    else{
      Lnewx = c(Lnewx,list(Lx[[i]]))
      Lnewt_x = c(Lnewt_x,list(Lt_x[[i]]))
    }
  }
  ty= y[isNewsub==0]
  if(is.na(optns_x)){
    optns_x = SetOptions(Lx_e,Lt_x_e,optns = list(methodSelectK='BIC'))
  }
  verbose_x = optns_x$verbose
  if(verbose_x == TRUE){
    print('Obtain functional object for x:\n')
  }
  xfpca = FPCA(Lx_e, Lt_x_e, optns_x) 
  xitrain = xfpca$xiEst
  
  zm = 50;
  ysort = sort(y)
  zgrid= seq(ysort[5], ysort[nn-5], length.out=zm)
  K_x = ncol(xitrain)
  b = matrix(0,zm,K_x+1)
  for(i in 1:zm){
    z=zgrid[i]
    ytrain =zeros(1,nn)
    ytrain[ty<=z]=1
    glmfit= glm(ytrain ~ xitrain, family = binomial(link= 'logit'))
    b[i,] = glmfit$coefficients
  }
  xi = matrix(0,n, K_x)
  xi[1:nn,] = xitrain
  if(length(newx) > 0){
    newpcx = FPCApred(xfpca,Lnewx,Lnewt_x)
    xi[(nn+1):n,] = newpcx
  }
  
  m = 100
  outF = seq(ysort[5], ysort[nn-5], length.out=m)
  predF = zeros(n, m)
  nQ = length(outQ)
  predQ = matrix(0,n, nQ)
  
  bw = 4*(zgrid[2]-zgrid[1])
  for(i in 1:n){
    CDF = 1:zm
    for(j in 1:zm){
      eta = b[j,] %*% xi[i,]
      CDF[j] = exp(eta)/(1+exp(eta))
    }
    win=array(1,zm);
    predF =lwls(bw,'epan', win, zgrid, CDF, outF)
  }
  
  for(i in 1:n){
    for(j in 1:nQ){
      if(predF[i, m] <= outQ[j]){
        predQ[i, j]= outF[m]
      }
      else{
        predQ[i, j]= outF[which(predF[i,] > outQ[j])[1]]
      }
    }
  }
  return(list(outQ,outF,b))
}