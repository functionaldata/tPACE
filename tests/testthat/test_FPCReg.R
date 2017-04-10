library(testthat)
#devtools::load_all('D:\\tPACE\\tPACE')

test_that ('Simple Dense Case works', {
set.seed(1000)
n <- 200 #number of subjects
ngrids <- 101 #number of grids in [0,1] for X(s)
ngridt <- 101 #number of grids in [0,1] for Y(t)
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
trueBetaX1Y <- sin(2 * pi * grids) %*% t(cos(2 * pi * gridt))
expect_equal(mean((estiBetaX1Y_Dense - trueBetaX1Y)^2), 0, tol = 1e-2)
})