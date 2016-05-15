devtools::load_all()
library(testthat)
##options(error=recover)

trueLam <- 4 / ((2 * (1:50) - 1 ) * pi) ^ 2

# set.seed(1)
# n <- 100
# pts <- seq(0, 1, by=0.05)
# samp3 <- Wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
# samp3 <- Sparsify(samp3, pts, 10)
# res <- FPCA(samp3$Ly, samp3$Lt, list(dataType='Sparse', useBins=TRUE))
# res$lambda / trueLam[1:length(res$lambda)]
# res$sigma2

# CreateCovPlot(res, 'Smoothed', FALSE)

test_that('error == FALSE works', {
  set.seed(1)
  n <- 100
  M <- 51
  pts <- seq(0, 1, length.out=M)
  mu <- rep(0, length(pts))
  sampDense <- Wiener(n, pts)
  samp4 <- Sparsify(sampDense, pts, 10)
  res4NE <- FPCA(samp4$Ly, samp4$Lt, list(error=FALSE, nRegGrid=M, FVEthreshold=1))
  res4Err <- FPCA(samp4$Ly, samp4$Lt, list(error=TRUE, nRegGrid=M, FVEthreshold=1))

  # CreatePathPlot(res4Err, 11:15)
  # CreatePathPlot(res4NE, 11:15)
  expect_equal(res4NE$sigma2, NULL)
  expect_true(mean((fitted(res4NE) - sampDense)^2) > 
              mean((fitted(res4Err) - sampDense)^2))
})

test_that('Truncation works for FPCA Wiener process', {
  set.seed(1)
  n <- 100
  pts <- seq(0, 1, by=0.01)
  mu <- rep(0, length(pts))
  samp4 <- Wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
  samp4 <- Sparsify(samp4, pts, 10)
  samp5 <- samp4
  samp4$Ly[[1]] <- samp4$Lt[[1]] <- c(0, 1)
  pTrunc <- SetOptions(samp4$Ly, samp4$Lt, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0.03, 0.97), verbose=TRUE))
  pNoTrunc <- SetOptions(samp4$Ly, samp4$Lt, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))
  set.seed(1); res4 <- FPCA(samp4$Ly, samp4$Lt, pTrunc)
  set.seed(1); res5NoTrunc <- FPCA(samp5$Ly, samp5$Lt, pNoTrunc)
  set.seed(1); res4NoTrunc <- FPCA(samp4$Ly, samp4$Lt, pNoTrunc)

  expect_equal(res4[c('sigma2', 'bwMu', 'bwCov')], res4NoTrunc[c('sigma2', 'bwMu', 'bwCov')])
  expect_equal(min(max(abs(res4$xiEst[-1, 1] - res4NoTrunc$xiEst[-1, 1])), max(abs(res4$xiEst[-1, 1] + res4NoTrunc$xiEst[-1, 1]))), 0, tol=0.5)
  expect_equal(min(max(abs(res5NoTrunc$xiEst[-1, 1] - res4NoTrunc$xiEst[-1, 1])), max(abs(res5NoTrunc$xiEst[-1, 1] + res4NoTrunc$xiEst[-1, 1]))), 0, tol=0.05)
  expect_equal(nrow(res4$xiEst), nrow(res4NoTrunc$xiEst))
  expect_equal(length(res4$xiVar), length(res4NoTrunc$xiVar))
})

test_that('Missing values work for FPCA Wiener process', {
  set.seed(1)
  n <- 200
  pts <- seq(0, 1, by=0.01)
  mu <- rep(0, length(pts))
  samp4 <- Wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
  samp4 <- Sparsify(samp4, pts, 10)
  samp4$Ly[[1]] <- samp4$Lt[[1]] <- c(0, 1)
  pNoTrunc <- SetOptions(samp4$Ly, samp4$Lt, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))

  samp4$Ly[[2]] <- samp4$Lt[[2]] <- c(0.1, 0.2, 0.5)
  set.seed(1); res4 <- FPCA(samp4$Ly, samp4$Lt, pNoTrunc)
  
  samp4$Ly[[2]] <- c(NA, 0.2, 0.5)
  set.seed(1); res4NaN <- FPCA(samp4$Ly, samp4$Lt, pNoTrunc)

  samp4$Ly[[2]] <-  c(0.1, 0.2, 0.5)
  samp4$Lt[[2]] <-  c(NA, 0.2, 0.5)
  set.seed(1); res4NaN2 <- FPCA(samp4$Ly, samp4$Lt, pNoTrunc)

  expect_equal(res4[c('sigma2', 'bwMu', 'bwCov')], res4NaN[c('sigma2', 'bwMu', 'bwCov')], tol=1e-7)
  expect_equal(nrow(res4$xiEst), nrow(res4NaN$xiEst))
  expect_equal(length(res4$xiVar), length(res4NaN$xiVar))
  expect_equal(res4NaN$inputData$Ly[[2]], c(0.2,0.5))
  expect_equal(   res4$inputData$Ly[[2]], c(0.1,0.2,0.5))
  expect_equal(res4NaN2[c('sigma2', 'bwMu', 'lambda')], res4NaN[c('sigma2', 'bwMu', 'lambda')], tol=1e-7)

})
test_that('User provided mu and cov for simple example',{

  set.seed(123)
  N = 200;   
  M = 100;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
 
  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) s  + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  ef <- matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)
  ev <- c(5, 2)^2
  covTrue <- ef %*% diag(ev) %*% t(ef)
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(sqrt(ev))
 
  # Create Y_true
  yTrue = Ksi %*% t(ef) + t(matrix(rep(meanFunct(s),N), nrow=M))
  
  # Create sparse sample  
  # Each subject has one to five readings (median: 3);
  ySparse = Sparsify(yTrue, s, c(2:5))
    
  # Give your sample a bit of noise 
  ySparse$yNoisy = lapply( ySparse$Ly, function(x) x +  1 * rnorm(length(x))) 
  
  # Do FPCA on this sparse sample
  FPCAsparseA = FPCA(ySparse$yNoisy,ySparse$Lt, optns = 
                    list(userMu = list(t=s,mu= meanFunct(s)), userCov = list(t=s,cov= covTrue) ))
  
  expect_equal(FPCAsparseA$sigma2, 1, tolerance=0.1)
  expect_equal( spline(y = FPCAsparseA$mu, x = FPCAsparseA$workGrid, xout = FPCAsparseA$obsGrid)$y, 
                expected = meanFunct(s), tolerance = 1e-3)
  expect_equal( FPCAsparseA$lambda, expected=ev, tolerance = 1e-1)
  expect_equal( abs(cor( FPCAsparseA$xiEst[,1], Ksi[,1])) > 0.9, TRUE)
})

test_that('User provided mu, cov, and sigma2',{

  set.seed(123)
  N = 200;   
  M = 90;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
 
  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) s  + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(c(5,2))
 
  # Create Y_true
  yTrue = Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))
  
  # Create sparse sample  
  # Each subject has one to five readings (median: 3);
  ySparse = Sparsify(yTrue, s, c(2:5))
    
  # Give your sample a bit of noise 
  ySparse$yNoisy = lapply( ySparse$Ly, function(x) x +  0.025*rnorm(length(x))) 
  
  userSigma2 <- 0.1
  # Do FPCA on this sparse sample
  FPCAsparseA = FPCA(ySparse$yNoisy,ySparse$Lt, optns = 
                    list(userMu = list(t=s,mu= meanFunct(s)), userSigma2=0.1 ))
    
  expect_equal( FPCAsparseA$sigma2, userSigma2)
  expect_equal( sqrt(FPCAsparseA$lambda[1:2]), expected=c(5,2), tolerance = 1e-1)
  expect_equal( abs(cor( FPCAsparseA$xiEst[,1], Ksi[,1])) > 0.95, TRUE)
  expect_equal( abs(cor( FPCAsparseA$xiEst[,2], Ksi[,2])) > 0.90, TRUE)
})

test_that('Case where one component should be returned',{

  set.seed(123)
  N = 111;   
  M = 81;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
 
  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) s  + 1 
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(c(5,2))
 
  # Create Y_true
  yTrue = Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))
 
  # Create sparse sample  
  # Each subject has one to five readings (median: 3);
  ySparse = Sparsify(yTrue, s, c(1:5))    
  FPCAsparseA = FPCA(ySparse$Ly,ySparse$Lt, optns = 
    list( FVEthreshold = 0.4, userMu = list(t=s,mu= meanFunct(s)) ) )

  expect_equal(spline(y = FPCAsparseA$mu, x = FPCAsparseA$workGrid, xout = FPCAsparseA$obsGrid)$y,
               expected = meanFunct(s), tolerance = 1e-3)
  expect_equal( length(FPCAsparseA$lambda), 1)

})

test_that('GetCovDense with noise, get sigma2', {
  set.seed(1)
  n <- 200
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  sampTrue <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE)
  samp <- sampTrue + rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- MakeFPCAInputs(tVec=pts, yVec=samp)

  resNoerr <- FPCA(tmp$Ly, tmp$Lt, list(error=FALSE, dataType='Dense', lean=TRUE))
  resErr <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', shrink=TRUE, lean=TRUE))
  resErr$optns

  expect_equal(resErr$sigma2, sigma2, tolerance=1e-1)
# shrinkage should be a bit better
  expect_true(mean((fitted(resNoerr) - sampTrue) ^ 2) > 
              mean((fitted(resErr) - sampTrue) ^ 2))
})

test_that('GetCovDense with noise, known cov, get sigma2', {
  set.seed(1)
  n <- 200
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  sampTrue <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE)
  samp <- sampTrue + rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- MakeFPCAInputs(tVec=pts, yVec=samp)

  resErr <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', shrink=TRUE, lean=TRUE, userCov=list(t=pts, cov=outer(pts, pts, pmin))))

  expect_equal(resErr$sigma2, sigma2, tolerance=1e-1)
})

test_that('noisy dense data, cross-sectional mu/cov, use IN/CE score', {
  set.seed(1)
  n <- 200
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  sampTrue <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE)
  samp <- sampTrue + rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- MakeFPCAInputs(tVec=pts, yVec=samp)

  resErrCE <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', lean=TRUE, methodXi='CE'))
  resErrCEKnownSigma2 <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', lean=TRUE, methodXi='CE', userSigma2=sigma2))
  resErrIN <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', lean=TRUE, methodXi='IN'))

  expect_equal(resErrCE$xiEst[, 1], resErrIN$xiEst[, 1], tolerance=1e-2)
  expect_equal(resErrCE$xiEst[, 2], resErrIN$xiEst[, 2], tolerance=5e-2)
  expect_equal(resErrCE$xiEst[, 3], resErrIN$xiEst[, 3], tolerance=1e-1)
  expect_equal(resErrCEKnownSigma2$xiEst[, 1], resErrIN$xiEst[, 1], tolerance=1e-2)
})

test_that('noisy dense data, smooth mu/cov, use IN/CE score', {
  set.seed(1)
  n <- 200
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  sampTrue <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE)
  samp <- sampTrue + rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- MakeFPCAInputs(tVec=pts, yVec=samp)

  resErrCSCE <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', lean=TRUE, methodMuCovEst='cross-sectional', methodXi='CE'))
  resErrCE <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', lean=TRUE, methodMuCovEst='smooth', methodXi='CE'))
  resErrIN <- FPCA(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense', lean=TRUE, methodMuCovEst='smooth', methodXi='IN'))

  trueCov <- outer(pts, pts, pmin)
  expect_true(max(abs(resErrCE$smoothedCov - trueCov)) < 0.1)
  expect_true(max(abs(resErrCSCE$smoothedCov - trueCov)) < 0.25)
  expect_equal(resErrCE$xiEst[, 1], resErrIN$xiEst[, 1], tolerance=1e-2)
  expect_equal(resErrCE$xiEst[, 2], resErrIN$xiEst[, 2], tolerance=6e-2)
  expect_equal(resErrCE$xiEst[, 3], resErrIN$xiEst[, 3], tolerance=1e-1)
})
