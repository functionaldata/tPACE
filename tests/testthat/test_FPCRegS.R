library(testthat)
library(MASS)

####################################
##Single predictor X1
set.seed(1000)
n = 300
N = 100
# eigenfunctions
phi11 <- function(t) sqrt(2)*sin(2*pi*t)
phi12 <- function(t) sqrt(2)*cos(2*pi*t)
sig <- diag(c(2,1.2))
scoreX <- mvrnorm(n,mu=rep(0,2),Sigma=sig)
scoreXTest <- mvrnorm(N,mu=rep(0,2),Sigma=sig)
# training set
gridX <- seq(0,1,length.out=21)
Lt <- Lx1 <- list()
# Lt <- Lx1 <- Lx2 <- list()
for (i in 1:n) {
  Lt[[i]] <- gridX
  Lx1[[i]] <- scoreX[i,1]*phi11(gridX) + scoreX[i,2]*phi12(gridX) + rnorm(1,0,0.01)
}
# test set
LtTest <- Lx1Test <- list()
for (i in 1:N) {
  LtTest[[i]] <- gridX
  Lx1Test[[i]] <- scoreXTest[i,1]*phi11(gridX) + scoreXTest[i,2]*phi12(gridX) + rnorm(1,0,0.01)
}
# dense observations
denseX1 <- matrix(nrow=n,ncol=length(gridX))
for (i in 1:n) {
  denseX1[i,] <- Lx1[[i]]
}
denseX1Test <- matrix(nrow=N,ncol=length(gridX))
for (i in 1:N) {
  denseX1Test[i,] <- Lx1Test[[i]]
}
denseX1Full <- rbind(denseX1,denseX1Test)
sig0 <- 0.2
betaCoef <- c(2,-1.5)
Y <- c(scoreX%*%betaCoef) + rnorm(n,0,sig0)
YTest <- c(scoreXTest%*%betaCoef) + rnorm(N,0,sig0)
YFull <- c(Y,YTest)

vars = list( X1 = list(Ly = append(Lx1,Lx1Test),Lt = append(Lt,LtTest) ), Y = Y)
isNewSub_DS = c(rep(0,length(Lx1)) ,rep(1,length(Lx1Test)))
result_DS = FPCRegS(vars,isNewSub = isNewSub_DS)

tmp1 = Sparsify(denseX1, gridX, 3:5)
tmp2 = Sparsify(denseX1Test, gridX, 3:5)
vars2 = list( X1 = list(Ly = append( tmp1$Ly,tmp2$Ly ),Lt = append(tmp1$Lt, tmp2$Lt) ), Y = Y)
isNewSub_SS = c(rep(0,length(tmp1$Ly)) ,rep(1,length(tmp2$Lt )))
resultSS = FPCRegS(vars2,isNewSub = isNewSub_SS,methodSelect = list(method = "Basis",NumberOfBasis = 2))

test_that ('Simple Dense and Sparse Case works', {

	expect_equal(mean((result_DS$predictions - YTest)^2)/var(YTest), 0, tol = 1e-2)
	expect_equal(mean((resultSS$predictions - YTest)^2)/var(YTest), 0, tol = 0.1)
}


####################################
##Multiple predictor X1 and X2
set.seed(1000)

n<-500
N<-200 
# eigenfunctions
phi11 <- function(t) sqrt(2)*sin(2*pi*t)
phi12 <- function(t) sqrt(2)*sin(4*pi*t)
phi21 <- function(t) sqrt(2)*cos(2*pi*t)
phi22 <- function(t) sqrt(2)*cos(4*pi*t)

###
### predictor generation (FPC scores from two functional predictors)
###

sig <- matrix(c(2.0, 0.0, 0.5, -.2,
                0.0, 1.2, -.2, 0.3,
                0.5, -.2, 1.7, 0.0,
                -.2, 0.3, 0.0, 1.0),
              nrow=4,ncol=4)
# FPC scores

library(MASS)
scoreX <- mvrnorm(n,mu=rep(0,4),Sigma=sig)
scoreXTest <- mvrnorm(N,mu=rep(0,4),Sigma=sig)

# training set
gridX <- seq(0,1,length.out=21)
Lt <- Lx1 <- Lx2 <- list()
for (i in 1:n) {
  Lt[[i]] <- gridX
  Lx1[[i]] <- scoreX[i,1]*phi11(gridX) + scoreX[i,2]*phi12(gridX) + rnorm(1,0,0.01)
  Lx2[[i]] <- scoreX[i,3]*phi21(gridX) + scoreX[i,4]*phi22(gridX) + rnorm(1,0,0.01)
}

# test set
LtTest <- Lx1Test <- Lx2Test <- list()
for (i in 1:N) {
  LtTest[[i]] <- gridX
  Lx1Test[[i]] <- scoreXTest[i,1]*phi11(gridX) + scoreXTest[i,2]*phi12(gridX) + rnorm(1,0,0.01)
  Lx2Test[[i]] <- scoreXTest[i,3]*phi21(gridX) + scoreXTest[i,4]*phi22(gridX) + rnorm(1,0,0.01)
}

# dense observations
denseX1 <- denseX2 <- matrix(nrow=n,ncol=length(gridX))
for (i in 1:n) {
  denseX1[i,] <- Lx1[[i]]
  denseX2[i,] <- Lx2[[i]]
}

denseX1Test <- denseX2Test <- matrix(nrow=N,ncol=length(gridX))
for (i in 1:N) {
  denseX1Test[i,] <- Lx1Test[[i]]
  denseX2Test[i,] <- Lx2Test[[i]]
}

denseX1Full <- rbind(denseX1,denseX1Test)
denseX2Full <- rbind(denseX2,denseX2Test)


###
### response generation (functional linear model)
###

betaCoef <- c(2,-1.5,1.8,1.3)
Y <- c(scoreX%*%betaCoef) + rnorm(n,0,0.5)
YTest <- c(scoreXTest%*%betaCoef) + rnorm(N,0,0.5)
YFull <- c(Y,YTest)

###
### linear regression
###

vars3 = list(X1=list(Ly = append(Lx1,Lx1Test),Lt = append(Lt,LtTest)),
             X2=list(Ly = append(Lx2,Lx2Test),Lt = append(Lt,LtTest)), 
             Y = Y)
isNewSub_DM = c(rep(0,length(Lx1)),rep(1,length(Lx1Test)))
result_DM <- FPCRegS(vars3,isNewSub = isNewSub_DM)

X1.SP = Sparsify(denseX1, gridX, 3:5)
X2.SP = Sparsify(denseX2, gridX, 2:7)
X1Test.SP = Sparsify(denseX1Test, gridX, 3:5)
X2Test.SP = Sparsify(denseX2Test, gridX, 2:7)

vars4 = list(X1=list(Ly = append(X1.SP$Ly,X1Test.SP$Ly),Lt = append(X1.SP$Lt,X1Test.SP$Lt)),
             X2=list(Ly = append(X2.SP$Ly,X2Test.SP$Ly),Lt = append(X2.SP$Lt,X2Test.SP$Lt)), 
             Y = Y)
isNewSub_SM = c(rep(0,length(X1.SP$Lt)),rep(1,length(X1Test.SP$Lt )))
#result6.1 <- FPCRegS(vars6,isNewSub = isNewSub.6,varsOptns = list(list(methodBwCov = 'GMeanAndGCV'),list(methodBwCov = 'GMeanAndGCV')))
result_SM <- FPCRegS(vars6,isNewSub = isNewSub_SM)




test_that ('Multiple Dense and Sparse Case works', {

	expect_equal(mean((result_DM$predictions - YTest)^2)/var(YTest), 0, tol = 0.05)
	expect_equal(mean((result_SM$predictions - YTest)^2)/var(YTest), 0, tol = 0.5)
}


test_that ('eigen decomposition of Cov works',{
expect_equal( abs(result_DM$Eigen[1] - diag(sig)[1])/abs(diag(sig)[1]), 0, tol = 0.1)
expect_equal( sum(abs(result_DM$Eigen[1:4] - diag(sig))/sum(diag(sig))), 0, tol = 0.3)
})







