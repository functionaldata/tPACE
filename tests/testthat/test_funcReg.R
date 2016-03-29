library(MASS)
library(devtools)
library(testthat)
load_all('../RPACE/tPACE')
source('func.R')
source('funcReg.R')


# Y(t) = \beta_0(t) + \beta_1(t) X_1(t) + \beta_2(t) X_2(t) + \beta_3 Z_3 + \beta_4 Z_4 + \epsilon
# X_1(t) = \mu_{X_1}(t) + Z_1 \phi_1(t); X_2(t) = \mu_{X_2}(t) + Z_2 \phi_2(t);
# \mu_{X_1}(t) = 2t, \mu_{X_2}(t) = -2t 
# (Z_1, \dots, Z_4) \sim N(0, \Sigma). \Sigma = 1/4 + diag(3/4)
# \phi_1(t) = 1, \phi_2(t) = 2.5 + cos(\pi t), t \in [0, 1]. So $cov(X_1(s), X_2(t)) = 1/4 \cos(\pi t).$
# \epsilon \sim N(0, \sigma^2).
# \beta_1 = \beta_2 = \beta_3 = \beta_4 = 1.

set.seed(1)
n <- 100
nGridIn <- 200
sparsity <- 5:10 # must have length > 1
bw <- 0.1
kern <- 'gauss'
T <- round(seq(0, 1, length.out=nGridIn), 4)

Sigma <- 1 / 4 + diag(3 / 4, 4)
mu <- T * 2
sigma <- 1

beta_0 <- 0
beta_1 <- 1
beta_2 <- 1
beta_3 <- 1
beta_4 <- 1

Z <- mvrnorm(n, rep(0, 4), Sigma)
X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(mu, n, nGridIn, byrow=TRUE)
X_2 <- Z[, 2, drop=FALSE] %*% matrix(2.5 + cos(pi * T), 1, nGridIn) - matrix(mu, n, nGridIn, byrow=TRUE)
epsilon <- rnorm(n, sd=sigma)
Y <- matrix(NA, n, nGridIn)
for (i in seq_len(n)) {
  Y[i, ] <- beta_0 + beta_1 * X_1[i, ] + beta_2 * X_2[i, ] + beta_3 * Z[i, 3] + beta_4 * Z[i, 4] + epsilon[i]
}
denseMuY <- colMeans(Y)
denseCovY <- cov(Y)

set.seed(1)
X_1sp <- Sparsify(X_1, T, sparsity)
set.seed(1)
X_2sp <- Sparsify(X_2, T, sparsity)
set.seed(1)
Ysp <- Sparsify(Y, T, sparsity)
Tout <- round(seq(min(T), 1, by=0.05), 2)
vars <- list(X_1=X_1sp, X_2=X_2sp, Z_3=Z[, 3], Z_4=Z[, 4], Y=Ysp)

test_that('Scaler-scaler cov', 
  expect_equal(as.numeric(uniCov(Z[, 3], Z[, 4])), cov(Z[, 3], Z[, 4]))
)

test_that('Scaler-function cov = Function-scaler cov', {
  expect_equal(uniCov(X_1sp, Z[, 3], bw, Tout), 
    t(uniCov(Z[, 3], X_1sp, bw, Tout)))
})
system.time(
cov12 <- uniCov(X_1sp, X_2sp, bw, Tout)
)
cov12_rd <- uniCov(X_1sp, X_2sp, bw, Tout, rmDiag=TRUE)
cov12_1D <- uniCov(X_1sp, X_2sp, bw, Tout, use1D=TRUE)
cov13 <- uniCov(X_1sp, Z[, 3], bw, Tout)
cov21 <- uniCov(X_2sp, X_1sp, bw, Tout)
cov11 <- uniCov(X_1sp, X_1sp, bw, Tout)
cov11_rd <- uniCov(X_1sp, X_1sp, bw, Tout, rmDiag=TRUE)
# diag(cov11)
# diag(cov11_rd)
# diag(cov11 - cov11_rd)
cov22 <- uniCov(X_2sp, X_2sp, bw, Tout)
cov22_rd <- uniCov(X_2sp, X_2sp, bw, Tout, rmDiag=TRUE)
# diag(cov22)
# diag(cov22_rd)
# diag(cov22 - cov22_rd)

test_that('Function-function cov works', {
  expect_equal(cov12, t(cov21))
  # rgl::persp3d(Tout, Tout, cov12, col='blue', xlab='X(s)', ylab='Y(t)')
  # rgl::persp3d(Tout, Tout, cov21, col='blue', xlab='X(s)', ylab='Y(t)')
# x-direction is close to constant.
  # expect_true(mean(apply(cov12, 2, sd)) < 0.1)

# y-direction is close to 1/4 (1.5 + cos(pi t))
  expect_true(sqrt(mean(colMeans(cov12) - 1 / 4 * (1.5 + cos(pi * Tout)))^2) < 0.2)
  
# 1D and 2D smoother is similar
  expect_equal(diag(cov12), diag(cov12_1D), tolerance=0.2)
})

covAll <- MvCov(vars, bw, Tout, kern)
covAllNoError <- MvCov(vars, bw, Tout, kern, measurementError=FALSE)

test_that('Multi-function/scaler cov works', {
# The cov(x, y) and cov(y, x) is symmetric.
  expect_equal(diag(3), 2 * diag(3), tolerance=Inf)
  expect_equal(as.numeric(cov12), as.numeric(covAll[, , 1, 2]))
  expect_equal(covAll[, , 3, 1], t(covAll[, , 1, 3]))
  expect_equal(covAll[, , 4, 3], t(covAll[, , 3, 4]))
  expect_equal(covAll[, , 1, 2], t(covAll[, , 2, 1]))
  expect_equal(covAll[, , 1, 2], covAllNoError[, , 1, 2])
  tmp <- max(abs(covAll[, , 1, 1] - covAllNoError[, , 1, 1]))
  expect_true(tmp > 0.01 && tmp < 0.35)

  # rgl::persp3d(Tout, Tout, covAll[, , 1, 1], col='blue', xlab='X(s)', ylab='Y(t)')
  # rgl::persp3d(Tout, Tout, covAll[, , 1, 1], col='blue', xlab='X(s)', ylab='Y(t)')
})

demeanRes <- demean(vars, bw, kern)
varsDemean <- demeanRes[['xList']]
muDemean <- demeanRes[['muList']]
covAllDemean <- MvCov(varsDemean, bw, Tout, kern, center=FALSE)

test_that('demean works', {
  expect_equal(demeanRes[['muList']][['Z_3']], mean(vars[['Z_3']]))
  expect_equal(demeanRes[['muList']][['Z_4']], mean(vars[['Z_4']]))
  expect_equal(uniCov(varsDemean[['X_1']], varsDemean[['X_2']], bw, Tout, center=FALSE), cov12)
  expect_equal(uniCov(varsDemean[['X_1']], varsDemean[['Z_3']], bw, Tout, center=FALSE), cov13)
  expect_equal(covAll, covAllDemean)
})

withError2D <- mvConReg(vars, bw, Tout)
withError1D <- mvConReg(vars, bw, Tout, diag1D='cross')
noError2D <- mvConReg(vars, bw, Tout, measurementError=FALSE)
noError1D <- mvConReg(vars, bw, Tout, measurementError=FALSE, diag1D='all')

# # Minimal eigenvalues sometimes smaller than 0.
# minLambda <- sapply(seq_along(Tout), function(i) {
  # min(eigen(noError1D[['cov']][i, i, 1:4, 1:4])[['values']])
# })
# minLambda2D <- sapply(seq_along(Tout), function(i) {
  # min(eigen(noError2D[['cov']][i, i, 1:4, 1:4])[['values']])
# })

test_that('1D and 2D covariance estimates are similar', {
  expect_equal(withError2D[['beta']], withError1D[['beta']], tolerance=0.2)
  expect_equal(noError2D[['beta']], noError1D[['beta']], tolerance=0.1)
})

withError2DRect <- mvConReg(vars, bw, Tout, kern='rect')
withError1DRect <- mvConReg(vars, bw, Tout, diag1D='cross', kern='rect')
noError2DRect <- mvConReg(vars, bw, Tout, measurementError=FALSE, kern='rect')
noError1DRect <- mvConReg(vars, bw, Tout, measurementError=FALSE, diag1D='all', kern='rect')
withError2DEpan <- mvConReg(vars, bw, Tout, kern='epan')
withError1DEpan <- mvConReg(vars, bw, Tout, diag1D='cross', kern='epan')
noError2DEpan <- mvConReg(vars, bw, Tout, measurementError=FALSE, kern='epan')
noError1DEpan <- mvConReg(vars, bw, Tout, measurementError=FALSE, diag1D='all', kern='epan')

test_that('Different kernel type works', {
# expect_equal(withError2DRect[['beta']], withError2D[['beta']], tolerance=0.01)
# expect_equal(withError1DRect[['beta']], withError1D[['beta']], tolerance=0.01)
# expect_equal(noError2DRect[['beta']], noError2D[['beta']], tolerance=0.01)
# expect_equal(noError1DRect[['beta']], noError1D[['beta']], tolerance=0.01)
  expect_equal(withError2DRect[['beta']], withError2DEpan[['beta']], tolerance=0.2)
  expect_equal(withError1DRect[['beta']], withError1DEpan[['beta']], tolerance=0.2)
  expect_equal(noError2DRect[['beta']], noError2DEpan[['beta']], tolerance=0.2)
  expect_equal(noError1DRect[['beta']], noError1DEpan[['beta']], tolerance=0.2)
})

fpcaX1 <- FPCA(X_1sp[['yList']], X_1sp[['tList']], list(userBwMu=bw, userBwCov=bw))
fpcaX2 <- FPCA(X_2sp[['yList']], X_2sp[['tList']], list(userBwMu=bw, userBwCov=bw))
fpcaY <- FPCA(Ysp[['yList']], Ysp[['tList']], list(userBwMu=bw, userBwCov=bw))
FPCAlist <- list(Y=fpcaY, X_1=fpcaX1, X_2=fpcaX2)
imputeRes <- imputeConReg(FPCAlist, Z[, 3:4], Tout)
test_that('imputation and 1D alpha beta estimates are similar', {
  expect_equal(imputeRes[['alpha']], noError1D[['alpha']], scale=1, tolerance=0.5)
  expect_equal(unname(imputeRes[['beta']]), noError1D[['beta']], scale=1, tolerance=0.5)
})

## subsetVars
test_that('subseting covariates is fine', {
  subVars <- subsetVars(vars, 1)
  subVarsTF <- subsetVars(vars, c(TRUE, rep(FALSE, n - 1)))
  expect_equal(subVars, subVarsTF)
  expect_equal(subVars[['X_1']][['tList']][[1]], X_1sp[['tList']][[1]])
  expect_equal(length(subVars[['X_1']][['tList']]), 1)
  expect_equal(length(subVars[['Z_3']]), 1)
})
