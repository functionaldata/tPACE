# devtools::load_all()
library(testthat)
# library(scatterplot3d)

try(silent=TRUE, load(system.file('testdata', 'InputFormMllwlskInCpp.RData', package='fdapace')))
IN <- InputFormMllwlskInCpp


test_that('Lwls2DDeriv original curve is correct', {
  A1 <- Lwls2D(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid, crosscov=TRUE)
  A2 <- Lwls2DDeriv(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid, crosscov=TRUE)
  B1 <- Lwls2D(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, crosscov=TRUE)
  B2 <- Lwls2DDeriv(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, crosscov=TRUE)
  C1 <- Lwls2D(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid, crosscov=TRUE)
  C2 <- Lwls2DDeriv(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid, crosscov=TRUE)
  D1 <- Lwls2D(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, crosscov=TRUE)
  D2 <- Lwls2DDeriv(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, crosscov=TRUE)
  expect_equal(A1, A2)
  expect_equal(B1, B2)
  expect_equal(C1, C2)
  expect_equal(D1, D2)
})


test_that('bilinear regression function', {
  n <- 500
  bw <- 0.7
  outGrid <- seq(-1, 1, by=0.05)
  sigma <- 0.1
  f <- function(x, y) -(2 * x + y)

  set.seed(1)
  xin <- matrix(runif(2 * n, -1, 1), ncol=2)
  yin <- apply(xin, 1, function(x) f(x[1], x[2]))
  inGrid <- seq(-1, 1, length.out=floor(sqrt(n)))
  xinReg <- expand.grid(inGrid, inGrid)
  yinReg <- apply(xinReg, 1, function(x) f(x[1], x[2]))

  # scatterplot3d(xin[, 1], xin[, 2], yin)
  # scatterplot3d(xinReg[, 1], xinReg[, 2], yinReg)

  # Noisy observations
  bwn <- 1
  yinNoisy <- yin + rnorm(length(yin), sd=sigma)
  yinRegNoisy <- yinReg + rnorm(length(yinReg), sd=sigma)
  # rgl::persp3d(inGrid, inGrid, yinRegNoisy, xlab='t1', ylab='t2')

  for (kern in c('gauss', 'epan')) {
    # Noiseless
    # Naming: results{npoly}{nder1}{nder2}
    resOld <- Lwls2D(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, crosscov=TRUE)
    d100 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=0, crosscov=TRUE)
    d110 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=1, nder2=0)
    d101 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=1)
    d211 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=2, nder1=1, nder2=1)
    # rgl::persp3d(outGrid, outGrid, d100)

    expect_equal(d100, sapply(outGrid, function(y) sapply(outGrid, function(x) f(x, y))))
    expect_equal(resOld, d100)
    expect_equal(d110, matrix(-2, nrow(d110), ncol(d110)))
    expect_equal(d101, matrix(-1, nrow(d101), ncol(d101)))
    expect_equal(d211, matrix(0, nrow(d211), ncol(d211)))

    # Noisy
    d100n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=0, crosscov=TRUE)
    d110n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=1, nder2=0)
    d101n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=1)
    # d111n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=1, nder2=1)
    d211n <- Lwls2DDeriv(1, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=2, nder1=1, nder2=1)
    # rgl::persp3d(outGrid, outGrid, d100n)
    # rgl::persp3d(outGrid, outGrid, d110n)
    # rgl::persp3d(outGrid, outGrid, d101n)
    # # rgl::persp3d(outGrid, outGrid, d111n)
    # rgl::persp3d(outGrid, outGrid, d211n)

    expect_equal(d100n, sapply(outGrid, function(y) sapply(outGrid, function(x) f(x, y))), tolerance=1e-2)
    expect_equal(d110n, matrix(-2, nrow(d110n), ncol(d110n)), tolerance=5e-2)
    expect_equal(d101n, matrix(-1, nrow(d101n), ncol(d101n)), tolerance=5e-2)
    expect_equal(d211n, matrix(0, nrow(d211n), ncol(d211n)), tolerance=1e-1)
  }
})


test_that('biquadratic regression function', {
  n <- 5000
  outGrid <- seq(-1, 1, by=0.05)
  sigma <- 1
  f <- function(x, y) -(2 * x + y)^2

  set.seed(1)
  xin <- matrix(runif(2 * n, -1, 1), ncol=2)
  yin <- apply(xin, 1, function(x) f(x[1], x[2]))
  inGrid <- seq(-1, 1, length.out=floor(sqrt(n)))
  xinReg <- expand.grid(inGrid, inGrid)
  yinReg <- apply(xinReg, 1, function(x) f(x[1], x[2]))

  # scatterplot3d(xin[, 1], xin[, 2], yin)
  # scatterplot3d(xinReg[, 1], xinReg[, 2], yinReg)

  # Noisy observations
  yinNoisy <- yin + rnorm(length(yin), sd=sigma)
  yinRegNoisy <- yinReg + rnorm(length(yinReg), sd=sigma)
  # rgl::persp3d(inGrid, inGrid, yinRegNoisy, xlab='t1', ylab='t2')

  for (kern in c('epan', 'gauss')) {
    if (kern == 'epan') {
      bw <- 0.2
    } else if (kern == 'gauss') {
      bw <- 0.1
    }
    bwn <- bw * 2 # bandwidth for noisy case
    
    # Noiseless
    # Naming: results{npoly}{nder1}{nder2}
    bw <- 0.1
    bwn <- 0.2
    resOld <- Lwls2D(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, crosscov=TRUE)
    d100 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=0, crosscov=TRUE)
    d110 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=1, nder2=0)
    d101 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=1)
    d211 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=2, nder1=1, nder2=1)
    # rgl::persp3d(outGrid, outGrid, d100)

    expect_equal(d100, 
      sapply(outGrid, function(y) sapply(outGrid, function(x) f(x, y))),
      tolerance=1e-1)
    expect_equal(resOld, d100)
    expect_equal(d110, 
      sapply(outGrid, function(y) sapply(outGrid, function(x) -4 * (2 * x + y))),
      tolerance=1e-1)
    expect_equal(d101, 
      sapply(outGrid, function(y) sapply(outGrid, function(x) -2 * (2 * x + y))),
      tolerance=1e-1)
    expect_equal(d211, matrix(-4, nrow(d211), ncol(d211)), tolerance=1e-5)


    d100n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=0, crosscov=TRUE)
    d110n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=1, nder2=0)
    d101n <- Lwls2DDeriv(bwn, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=1)
    d211n <- Lwls2DDeriv(1, kern, xin, yinNoisy, xout1=outGrid, xout2=outGrid, npoly=2, nder1=1, nder2=1)
    # rgl::persp3d(outGrid, outGrid, d100n)
    # rgl::persp3d(outGrid, outGrid, d110n)
    # rgl::persp3d(outGrid, outGrid, d101n)
    # rgl::persp3d(outGrid, outGrid, d211n)

    expect_equal(d100n, 
      sapply(outGrid, function(y) sapply(outGrid, function(x) f(x, y))),
      tolerance=1e-1)
    expect_equal(d110n, 
      sapply(outGrid, function(y) sapply(outGrid, function(x) -4 * (2 * x + y))),
      tolerance=5e-1)
    expect_equal(d101n, 
      sapply(outGrid, function(y) sapply(outGrid, function(x) -2 * (2 * x + y))),
      tolerance=5e-1)
    expect_equal(d211n, matrix(-4, nrow(d211n), ncol(d211n)), tolerance=1e-1)
  }
})


test_that('biexponential regression function', {
  sqrtn <- 50
  n <- sqrtn^2
  inGrid <- outGrid <- seq(-1, 1, length.out=sqrtn)
  f <- function(x, y) exp(-2 * x + y)

  set.seed(1)
  xin <- as.matrix(expand.grid(inGrid, inGrid))
  yin <- f(xin[, 1], xin[, 2])
  inGrid <- seq(-1, 1, length.out=floor(sqrt(n)))
  # persp3d(inGrid, inGrid, yin, col='white', xlab='x', ylab='y')

  kern <- 'epan'
  bw <- 0.1

  # Noiseless
  # Naming: results{npoly}{nder1}{nder2}
  resOld <- Lwls2D(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, crosscov=TRUE)
  d100 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=0, crosscov=TRUE)
  d110 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=1, nder2=0)
  d101 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=1, nder1=0, nder2=1)
  d211 <- Lwls2DDeriv(bw, kern, xin, yin, xout1=outGrid, xout2=outGrid, npoly=2, nder1=1, nder2=1)

  expect_equal(d100, 
               matrix(yin, sqrtn, sqrtn),
               tolerance=1e-2)
  expect_equal(resOld, d100)
  expect_equal(d110, 
               sapply(outGrid, function(y) 
                      sapply(outGrid, function(x) -2 * exp(-2 * x + y))), 
               tolerance=1e-1)
  expect_equal(d101, 
               sapply(outGrid, function(y) 
                      sapply(outGrid, function(x) exp(-2 * x + y))),
               tolerance=1e-1)
  expect_equal(d211, sapply(outGrid, function(y) 
                            sapply(outGrid, function(x) -2 * exp(-2 * x + y))), tolerance=1e-1)

})
