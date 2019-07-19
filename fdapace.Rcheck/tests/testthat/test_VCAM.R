library(testthat)
#devtools::load_all()
library(MASS)
require(copula)||install.packages(copula)

MISE<-function(est,beta,timepoint){
  res<-0
  for(i in 2:length(timepoint)){
    res = res + (est[i]-beta[i])^2*(timepoint[i] - timepoint[i-1])
  }
  res
}
test_that('equidistance time point example works.',{
  set.seed(100)
  n <- 100
  d <- 2
  Lt <- list()
  Ly <- list()
  myCop <- normalCopula(param=0.5, dim = 2)
  myMvd <- mvdc(copula=myCop, margins=c("unif", "unif"),paramMargins=list(list(min = 0 , max = 1),list(min = 0 , max = 1) ))
  X <- rMvdc(n,myMvd)
  beta0 <- function(t) 1.5*sin(3*pi*(t+0.5))+4*t^3
  beta1 <- function(t) 3*(1-t)^2
  beta2 <- function(t) 4*t^3
  
  phi1 <- function(x) sin(2*pi*x)
  phi2 <- function(x) 4*x^3-1
  
  for (i in 1:n) {
    Lt[[i]] <- seq(0, 1, length.out = 21)
    Ly[[i]] <- beta0(Lt[[i]]) + beta1(Lt[[i]])*phi1(X[i,1]) + beta2(Lt[[i]])*phi2(X[i,2]) + rnorm(21,0,0.1)
    A = mvrnorm(1,rep(0,4),diag(c(1/4,1/9,1/16,1/25)))
    U = A[1]*sqrt(2)*cos(2*pi*Lt[[i]])+A[2]*sqrt(2)*sin(2*pi*Lt[[i]])+A[3]*sqrt(2)*cos(4*pi*Lt[[i]])+A[4]*sqrt(2)*sin(4*pi*Lt[[i]])
    Ly[[i]] <- Ly[[i]]+U
  }
  vcam <- VCAM(Lt,Ly,X)
  expect_lt( MISE(vcam$beta0Est,beta0(vcam$gridT),vcam$gridT), 0.05)
  expect_lt( MISE(vcam$betaEst[,1],beta1(vcam$gridT),vcam$gridT), 0.05)
  expect_lt( MISE(vcam$betaEst[,2],beta2(vcam$gridT),vcam$gridT), 0.05)
})

test_that('multiple covariates example works.',{
  set.seed(1071)
  n <- 100
  d <- 3
  Lt <- list()
  Ly <- list()
  myCop <- normalCopula(param=c(0.8,0.2,-0.4), dim = 3,'dispstr' = "un")
  myMvd <- mvdc(copula=myCop, margins=c("unif", "unif", "unif"),paramMargins=list(list(min = 0 , max = 1),list(min = 0 , max = 1),list(min = 0 , max = 1) ))
  X <- rMvdc(n,myMvd)
  beta0 <- function(t) 1.5*sin(3*pi*(t+0.5))+4*t^3
  beta1 <- function(t) 3*(1-t)^2
  beta2 <- function(t) 4*t^3
  beta3 <- function(t) sqrt(2)*cos(4*pi*t)+1
  phi1 <- function(x) sin(2*pi*x)
  phi2 <- function(x) 4*x^3-1
  phi3 <- function(x) sin(4*pi*x)
  for (i in 1:n) {
    Lt[[i]] <- seq(0, 1, length.out = 21)
    Ly[[i]] <- beta0(Lt[[i]]) + beta1(Lt[[i]])*phi1(X[i,1]) + beta2(Lt[[i]])*phi2(X[i,2]) + beta3(Lt[[i]])*phi3(X[i,3]) + rnorm(21,0,0.1)
    A = mvrnorm(1,rep(0,4),diag(c(1/4,1/9,1/16,1/25)))
    U = A[1]*sqrt(2)*cos(2*pi*Lt[[i]])+A[2]*sqrt(2)*sin(2*pi*Lt[[i]])+A[3]*sqrt(2)*cos(4*pi*Lt[[i]])+A[4]*sqrt(2)*sin(4*pi*Lt[[i]])
    Ly[[i]] <- Ly[[i]]+U
  }
  vcam <- VCAM(Lt,Ly,X)
  expect_lt(MISE(vcam$beta0Est,beta0(vcam$gridT),vcam$gridT) , 0.05)
  expect_lt(MISE(vcam$betaEst[,1],beta1(vcam$gridT),vcam$gridT), 0.05)
  expect_lt(MISE(vcam$betaEst[,2],beta2(vcam$gridT),vcam$gridT), 0.05)
  expect_lt(MISE(vcam$betaEst[,3],beta2(vcam$gridT),vcam$gridT) , 3)
})


test_that('Irregular design example works.',{
  set.seed(123)
  n <- 100
  d <- 2
  Lt <- list()
  Ly <- list()
  
  myCop <- normalCopula(param=0.5, dim = 2)
  myMvd <- mvdc(copula=myCop, margins=c("unif", "unif"),paramMargins=list(list(min = 0 , max = 1),list(min = 0 , max = 1) ))
  X <- rMvdc(n,myMvd)
  
  beta0 <- function(t) 1.5*sin(3*pi*(t+0.5))+4*t^3
  beta1 <- function(t) 3*(1-t)^2
  beta2 <- function(t) 4*t^3
  
  phi1 <- function(x) sin(2*pi*x)
  phi2 <- function(x) 4*x^3-1
  
  for (i in 1:n) {
    Ni <- sample(10:20,1)
    Lt[[i]] <- sort(runif(Ni,0,1))
    Ly[[i]] <- beta0(Lt[[i]]) + beta1(Lt[[i]])*phi1(X[i,1]) + beta2(Lt[[i]])*phi2(X[i,2]) + rnorm(Ni,0,0.1)
    A = mvrnorm(1,rep(0,4),diag(c(1/4,1/9,1/16,1/25)))
    U = A[1]*sqrt(2)*cos(2*pi*Lt[[i]])+A[2]*sqrt(2)*sin(2*pi*Lt[[i]])+A[3]*sqrt(2)*cos(4*pi*Lt[[i]])+A[4]*sqrt(2)*sin(4*pi*Lt[[i]])
    Ly[[i]] <- Ly[[i]]+U
  }
  vcam <- VCAM(Lt,Ly,X)
  expect_lt(MISE(vcam$beta0Est,beta0(vcam$gridT),vcam$gridT) , 0.1)
  expect_lt(MISE(vcam$betaEst[,1],beta1(vcam$gridT),vcam$gridT), 0.1)
  expect_lt(MISE(vcam$betaEst[,2],beta2(vcam$gridT),vcam$gridT), 0.1)
})
