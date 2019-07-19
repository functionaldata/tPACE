# devtools::load_all()

test_that('GetNormalisedSample output a homoscadestic sample', {
  set.seed(1)
  n <- 100
  M <- 51
  pts <- seq(0, 1, length.out=M)
  mu <- rep(0, length(pts))
  sampDense <- MakeGPFunctionalData(n, M, mu, K=1, basisType='sin', sigma=0.01)
  samp4 <- MakeFPCAInputs(tVec=sampDense$pts, yVec=sampDense$Yn)
  res4E <- FPCA(samp4$Ly, samp4$Lt, list(error=TRUE))
  sampN <- GetNormalisedSample(res4E, errorSigma=TRUE)
  sampN0 <- GetNormalisedSample(res4E, errorSigma=FALSE)

  # Cross-sectional standard deviation
  sdCr <- apply(simplify2array(sampN$Ly), 1, sd)
  sdCr0 <- apply(simplify2array(sampN0$Ly), 1, sd)
  expect_equal(sdCr[-c(1:2, (M-1):M)], rep(1, M - 4), tolerance=1e-4)
  expect_equal(sdCr0[-c(1:2, (M-1):M)], rep(1, M - 4), tolerance=1e-3)

  # CreatePathPlot(subset=1:20, inputData=samp4, obsOnly=TRUE, showObs=FALSE)
  # CreatePathPlot(subset=1:20, inputData=sampN, obsOnly=TRUE, showObs=FALSE)
  # CreatePathPlot(subset=1:20, inputData=sampN0, obsOnly=TRUE, showObs=FALSE)
})
