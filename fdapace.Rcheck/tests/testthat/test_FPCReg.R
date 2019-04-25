library(testthat)
#devtools::load_all('D:\\tPACE\\tPACE')

test_that ('Simple Dense and Sparse Case works', {
set.seed(1000)
n <- 300 #number of subjects
ngrids <- 51 #number of grids in [0,1] for X(s)
ngridt <- 51 #number of grids in [0,1] for Y(t)
grids <- seq(0, 1, length.out = ngrids) #regular grids in [0,1] for X(s)
gridt <- seq(0, 1, length.out = ngridt) #regular grids in [0,1] for Y(t)
#generate X
#{1, sqrt(2)*sin(2*pi*s), sqrt(2)*cos(2*pi*t)} are used to generate X.
eigenFun <- list(function(s){1 + 0 * s}, function(s){sqrt(2) * sin(2 * pi * s)}, function(s){sqrt(2) * cos(2 * pi * s)})
basisX <- sapply(eigenFun,function(x){x(grids)})
scoreX <- matrix(rnorm(n * 3), n, 3) #eigenvalues are assumed to be 1.
latentX <- scoreX %*% t(basisX)
measErrX <- sqrt(0.01) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
denseX <- latentX + measErrX
#generate Y
#beta(s,t) <- sin(2*pi*s)*cos(2*pi*t)
betaEigen1 <- function(t){f <- function(s){sin(2 * pi * s) * cos(2 * pi * t) * (1 + 0 * s)}; return(f)}
betaEigen2 <- function(t){f <- function(s){sin(2 * pi * s) * cos(2 * pi * t) * (sqrt(2) * sin(2 * pi * s))}; return(f)}
betaEigen3 <- function(t){f <- function(s){sin(2 * pi * s) * cos(2 * pi * t) * (sqrt(2) * cos(2 * pi * s))}; return(f)}
betaEigen <- list(betaEigen1, betaEigen2, betaEigen3) 
basisY <- array(0, c(ngridt, 3))
for (i in 1:3) {
	intbetaEigen <- function (t) {integrate(betaEigen[[i]](t), lower=0, upper=1)$value}
	basisY[,i] <- sapply(1:ngridt, function(x){intbetaEigen(gridt[x])})
	}
latentY <- scoreX %*% t(basisY)
measErrY <- sqrt(0.01) * matrix(rnorm(n * ngridt), n, ngridt) #0.05 is sigma^2
denseY <- latentY + measErrY

#Dense data
timeX <- t(matrix(rep(grids, n), length(grids), n))
timeY <- t(matrix(rep(gridt, n), length(gridt), n))
denseVars <- list(X1 = list(Ly = denseX, Lt = timeX), Y = list(Ly = denseY, Lt = timeY))
resuDense <- FPCReg(denseVars)
estiBetaX1Y_Dense <- resuDense$estiBeta$betaX1Y

#Sparse data with fixed bw to avoid sensitive issue of possible revised of FPCA() bw in near future
sparseX <- Sparsify(denseX,grids,1:10)
sparseY <- Sparsify(denseY,gridt,1:10)
sparseVars <- list(X1 = sparseX, Y = sparseY)
resuSparse <- FPCReg(sparseVars,varsOptns = list(X1=list(userBwMu = 0.05, userBwCov = 0.1),Y=list(userBwMu = 0.05, userBwCov = 0.1)),Kx=3)
estiBetaX1Y_Sparse <- resuSparse $estiBeta$betaX1Y

trueBetaX1Y <- sin(2 * pi * grids) %*% t(cos(2 * pi * gridt))

#Dense ratio of err and trun under L^2
expect_equal(mean((estiBetaX1Y_Dense - trueBetaX1Y)^2)/mean(trueBetaX1Y^2), 0, tol = 1e-2)
#Sparse ratio of err and trun under L^2
expect_equal(mean((estiBetaX1Y_Sparse - trueBetaX1Y)^2)/mean(trueBetaX1Y^2), 0, tol = 0.3)
})



#test for two predictors X1 and X2 
set.seed(1000)
#Model: E(Y(t)|X) = int(beta(s,t)*X(s))
n <- 300 #number of subjects
ngrids <- 101 #number of grids in [0,1] for X(s)
ngridt <- 51 #number of grids in [0,1] for Y(t)
grids <- seq(0, 1, length.out=ngrids) #regular grids in [0,1] for X(s)
gridt <- seq(0, 1, length.out=ngridt) #regular grids in [0,1] for Y(t)
 
#generate X
#{1, sqrt(2)*sin(2*pi*s), sqrt(2)*cos(2*pi*t)} are used to generate X.
eigenFun <- list(function(s){1 + 0 * s},function(s){sqrt(2) * sin(2*pi*s)},function(s){sqrt(2) * cos(2*pi*s)})
 
sig <- matrix(c(1.5, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.2, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.5, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
                nrow=6,ncol=6)

scoreX <- MASS::mvrnorm(n,mu=rep(0,6),Sigma=sig)
scoreX1 <- scoreX[,1:3]
scoreX2 <- scoreX[,4:6]

basisX1 <- sapply(eigenFun,function(x){x(grids)})
latentX1 <- scoreX1 %*% t(basisX1)
measErrX1 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
denseX1 <- latentX1 + measErrX1
 
basisX2 <- sapply(eigenFun,function(x){x(grids)})
latentX2 <- scoreX2 %*% t(basisX2)
measErrX2 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
denseX2 <- latentX2 + measErrX2

#generate Y
#beta(s, t) <- sin(2 * pi * s)*cos(2 * pi * t)
betaEigen1 <- function(t){f <- function(s){sin(2*pi*s) * cos(2*pi*t) * (1+0*s)};return(f)}
betaEigen2 <- function(t){f <- function(s){sin(2*pi*s) * cos(2*pi*t) * (sqrt(2)*sin(2*pi*s))}; return(f)}
betaEigen3 <- function(t){f <- function(s){sin(2*pi*s) * cos(2*pi*t) * (sqrt(2)*cos(2*pi*s))}; return(f)}
betaEigen <- list(betaEigen1, betaEigen2, betaEigen3) 
basisY <- array(0,c(ngridt, 3))
for(i in 1:3){
 	intbetaEigen <- function (t) {integrate(betaEigen[[i]](t), lower = 0, upper = 1)$value}
 	basisY[, i] <- sapply(1:ngridt, function(x){intbetaEigen(gridt[x])})
	}
#Note that in the case BetaX1Y = beta and BetaX2Y = -beta  
latentY <- scoreX1 %*% t(basisY) - scoreX2 %*% t(basisY)
measErrY <- sqrt(0.01) * matrix(rnorm(n*ngridt), n, ngridt) #0.01 is sigma^2
denseY <- latentY + measErrY

#======Dense data===============================================
timeX <- t(matrix(rep(grids, n),length(grids), n))
timeY <- t(matrix(rep(gridt, n),length(gridt), n))
denseVars <- list(X1 = list(Ly = denseX1, Lt = timeX),
                  X2 = list(Ly = denseX2, Lt = timeX),
                  Y=list(Ly = denseY,Lt = timeY))

p <- length(denseVars)-1
optns <- list(dataType = "Dense" ,error = 1, kernel='gauss', nRegGrid=51, useBinnedData='OFF')
brkX <- c(0, cumsum(c(ngrids,ngrids)))
dxMatrix <- dx(p, c(1,1), c(ngrids,ngrids), brkX)
gridNum <- c(ngrids,ngrids)
resdenseCov <- denseCov(p, denseVars, brkX, dxMatrix, gridNum, varsOptns= list(X1=optns,X2=optns,Y=optns))

test_that ('eigen decomposition of denseCov() works',{
expect_equal( abs(resdenseCov$eigenValue[1]-max(diag(sig)))/max(diag(sig)), 0, tol = 1e-1)
expect_equal( abs(sum(resdenseCov$eigenValue[1:6])-sum(diag(sig)))/sum(diag(sig)), 0, tol = 1e-1)
})

resuDense <- FPCReg(denseVars, method="FVE") 
estiBetaX1Y_Dense <- resuDense$estiBeta$betaX1Y
estiBetaX2Y_Dense <- resuDense$estiBeta$betaX2Y

trueBetaX1Y <- sin(2 * pi * grids) %*% t(cos(2 * pi * gridt))
test_that ('beta works in Dense with two predictors',{
expect_equal(mean((estiBetaX1Y_Dense - trueBetaX1Y)^2)/mean(trueBetaX1Y^2), 0, tol = 1e-2)
expect_equal(mean((estiBetaX2Y_Dense + trueBetaX1Y)^2)/mean(trueBetaX1Y^2), 0, tol = 1e-2)
})
 
#======Sparse data===============================================
sparsity = 5:10
sparseX1 <- Sparsify(denseX1, grids, sparsity)
sparseX2 <- Sparsify(denseX2, grids, sparsity)
sparseY <- Sparsify(denseY, gridt, sparsity)
sparseVars <- list(X1 = sparseX1, X2 = sparseX2, Y = sparseY)

ngrids <- 51
grid <- seq(0, 1, length.out = ngrids)
p <- length(denseVars)-1
brkX <- c(0, cumsum(c(ngrids,ngrids)))
dxMatrix <- dx(p, c(1,1), c(ngrids,ngrids), brkX)
gridNum <- c(ngrids,ngrids)
optns <- list(dataType = "Sparse" ,error = TRUE, kernel='gauss', nRegGrid=ngrids , useBinnedData='OFF')
demeanedRes <- demeanFuc(p, sparseVars, kern='gauss', varsOptns= list(X1=optns,X2=optns,Y=optns)) #Centered predictors. Using gauss for demeanFuc, but may be relaxed.
varsTrain <- demeanedRes[['xList']]
muList <- demeanedRes[['muList']]
ressparseCov <- sparseCov(p, sparseVars, brkX, varsOptns= list(X1=optns,X2=optns,Y=optns), sparseVars, muList, gridNum, dxMatrix, list(grid,grid))
 
test_that ('eigen decomposition of sparseCov() works',{
expect_equal( abs(ressparseCov$eigenValue[1]-max(diag(sig)))/max(diag(sig)), 0, tol = 0.3)
expect_equal( abs(sum(ressparseCov$eigenValue[1:6])-sum(diag(sig)))/sum(diag(sig)), 0, tol = 0.3)
})

resuSparse <- FPCReg(sparseVars, varsOptns = list(X1=list(userBwMu = 0.05, userBwCov = 0.1), X2=list(userBwMu = 0.05, userBwCov = 0.1), Y=list(userBwMu = 0.05, userBwCov = 0.1)) , Kx=6) 
estiBetaX1Y_Sparse <- resuSparse$estiBeta$betaX1Y 
estiBetaX2Y_Sparse <- resuSparse$estiBeta$betaX2Y

trueBetaX1Y <- sin(2 * pi * grid) %*% t(cos(2 * pi * grid))
test_that ('beta estimates works in sparse with two predictors',{
expect_equal(mean((estiBetaX1Y_Sparse - trueBetaX1Y)^2)/mean(trueBetaX1Y^2), 0, tol = 0.3)
expect_equal(mean((estiBetaX2Y_Sparse + trueBetaX1Y)^2)/mean(trueBetaX1Y^2), 0, tol = 0.3)
})