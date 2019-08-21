#setwd('/Users/kyunghee/Desktop/SBF')
#SBF_scripts <- list.files(pattern="*.R")
#for(i in 1:length(SBF_scripts)){
#  source(SBF_scripts[i])
#}

#library(Rcpp)
#setwd('/Users/kyunghee/Desktop/tPACE/src')
#sourceCpp('trapzRcpp.cpp')

library(testthat)

test_that(
  '(1) algorithm convergence including dimension checking, (2) theoretical estimation precision of true component functions.',
  {
    set.seed(100)
    n <- 100
    d <- 2
    X <- pnorm(matrix(rnorm(n*d),nrow=n,ncol=d)%*%matrix(c(1,0.6,0.6,1),nrow=2,ncol=2))
    
    f1 <- function(t) 2*(t-0.5)
    f2 <- function(t) sin(2*pi*t)
    
    Y <- f1(X[,1])+f2(X[,2])+rnorm(n,0,0.1)
    
    N <- 101
    x <- matrix(rep(seq(0,1,length.out=N),d),nrow=N,ncol=d)
    h <- c(0.12,0.08)
    
    sbfResult <- SBFitting(Y,x,X,h)
    
    fFit <- sbfResult$SBFit
    
    iterErr <- sbfResult$iterErr
    iterErrDiff <- sbfResult$iterErrDiff
    iterNum <- sbfResult$iterNum
    critErr <- sbfResult$critErr
    critErrDiff <- sbfResult$critErrDiff
    critNum <- sbfResult$critNum
    
    expect_true(sum(dim(fFit)!=c(N,d))==0)
    expect_true((iterErr<critErr)+(iterErrDiff<critErrDiff)+(iterNum<critNum)<3)
    
    estErr <- c(max(abs(f1(x[,1])-fFit[,1])),max(abs(f2(x[,2])-fFit[,2])))
    
    crit1 <- sqrt(log(n)/n/h^2)
    crit2 <- h
    
    crit <- apply(cbind(crit1,crit2),1,'max')
    
    expect_true(estErr[1]/crit[1]<1)
    expect_true(estErr[2]/crit[2]<1)
  }
)


test_that(
  '(1) algorithm convergence including dimension checking, (2) better prediction performance than constant or linear models.',
  {
    set.seed(100)
    n <- 100
    d <- 2
    X <- pnorm(matrix(rnorm(n*d),nrow=n,ncol=d)%*%matrix(c(1,0.6,0.6,1),nrow=2,ncol=2))
    
    f1 <- function(t) 2*(t-0.5)
    f2 <- function(t) sin(2*pi*t)
    
    Y <- f1(X[,1])+f2(X[,2])+rnorm(n,0,0.1)
    
    N <- n
    x <- X
    h <- c(0.12,0.08)
    
    sbfResult <- SBFitting(Y,x,X,h)
    
    fFit <- sbfResult$SBFit
    
    iterErr <- sbfResult$iterErr
    iterErrDiff <- sbfResult$iterErrDiff
    iterNum <- sbfResult$iterNum
    critErr <- sbfResult$critErr
    critErrDiff <- sbfResult$critErrDiff
    critNum <- sbfResult$critNum
    
    expect_true(sum(dim(fFit)!=c(N,d))==0)
    expect_true((iterErr<critErr)+(iterErrDiff<critErrDiff)+(iterNum<critNum)<3)
    
    mY <- sbfResult$mY
    
    mseSBF <- mean((Y - mY - apply(fFit,1,'sum'))^2)
    
    mseConst <- var(Y)
    mseLM <- mean(lm(Y~X)$resid^2)
    
    expect_true(mseSBF<mseConst)
    expect_true(mseSBF<mseLM)
  }
)