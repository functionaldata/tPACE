setwd('/Users/kyunghee/Desktop/SBF')
SBF_scripts<-list.files(pattern="*.R")
for(i in 1:length(SBF_scripts)){
  source(SBF_scripts[i])
}

library(Rcpp)
setwd('/Users/kyunghee/Desktop/tPACE/src')
sourceCpp('trapzRcpp.cpp')

library(testthat)

test_that(
  'estimation error = O(h) + O_P(sqrt(log(n)/(nh^2)))',
  {
    #set.seed(100)
    n<-100
    d<-2
    X<-pnorm(matrix(rnorm(n*d),nrow=n,ncol=d)%*%matrix(c(1,0.6,0.6,1),nrow=2,ncol=2))
    
    f1<-function(t) 2*(t-0.5)
    f2<-function(t) sin(2*pi*t)
    
    Y<-f1(X[,1])+f2(X[,2])+rnorm(n,0,0.1)
    
    N<-101
    x<-matrix(rep(seq(0,1,length.out=N),d),nrow=N,ncol=d)
    h<-c(0.12,0.08)
    
    SBF_result<-SBFitting(Y,x,X,h)
    f_fit<-SBF_result$SBFit
    
    est_err<-c(max(abs(f1(x[,1])-f_fit[,1])),max(abs(f2(x[,2])-f_fit[,2])))
    
    crit1<-sqrt(log(n)/n/h^2)
    crit2<-h
    
    crit<-apply(cbind(crit1,crit2),1,'max')
    
    expect_true(est_err[1]/crit[1]<0.25)
    expect_true(est_err[2]/crit[2]<0.25)
  }
)
