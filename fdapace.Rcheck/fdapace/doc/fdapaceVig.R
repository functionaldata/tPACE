## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  FPCAobj <- FPCA(Ly=yList, Lt=tList)

## ----eval=TRUE, echo=TRUE------------------------------------------------
  library(fdapace)
 
  # Set the number of subjects (N) and the
  # number of measurements per subjects (M) 
  N <- 200;
  M <- 100;
  set.seed(123)

  # Define the continuum
  s <- seq(0,10,length.out = M)

  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) s + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)

  # Create FPC scores
  Ksi <- matrix(rnorm(N*2), ncol=2);
  Ksi <- apply(Ksi, 2, scale)
  Ksi <- Ksi %*% diag(c(5,2))

  # Create Y_true
  yTrue <- Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

## ----eval=TRUE, echo=TRUE------------------------------------------------
  L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(s,N), t(yTrue))
  FPCAdense <- FPCA(L3$Ly, L3$Lt)

  # Plot the FPCA object
  plot(FPCAdense)

  # Find the standard deviation associated with each component
  sqrt(FPCAdense$lambda)
  

## ----eval=TRUE, echo=TRUE------------------------------------------------
  # Create sparse sample  
  # Each subject has one to five readings (median: 3)
  set.seed(123)
  ySparse <- Sparsify(yTrue, s, sparsity = c(1:5))

  # Give your sample a bit of noise 
  ySparse$yNoisy <- lapply(ySparse$Ly, function(x) x + 0.5*rnorm(length(x)))

  # Do FPCA on this sparse sample
  # Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
  # Smoothing is the main computational cost behind sparse FPCA
  FPCAsparse <- FPCA(ySparse$yNoisy, ySparse$Lt, list(plot = TRUE))


## ----eval=TRUE, echo=TRUE------------------------------------------------
 FPCAsparseMuBW5 <- FPCA(ySparse$yNoisy, ySparse$Lt, optns= list(userBwMu = 5))

## ----eval=TRUE, echo=TRUE------------------------------------------------
par(mfrow=c(1,2))
CreatePathPlot( FPCAsparse, subset = 1:3, main = "GCV bandwidth", pch = 16)
CreatePathPlot( FPCAsparseMuBW5, subset = 1:3, main = "User-defined bandwidth", pch = 16)

## ----eval=TRUE, echo=TRUE------------------------------------------------
 FPCAsparseRect <- FPCA(ySparse$yNoisy, ySparse$Lt, optns = list(kernel = 'rect')) # Use rectangular kernel

## ----eval=TRUE, echo=TRUE------------------------------------------------
SelectK( FPCAsparse, criterion = 'FVE', FVEthreshold = 0.95) # K = 2
SelectK( FPCAsparse, criterion = 'AIC') # K = 2

## ----eval=TRUE, echo=TRUE------------------------------------------------
fittedCurvesP0 <- fitted(FPCAsparse) # equivalent: fitted(FPCAsparse, derOptns=list(p = 0));
# Get first order derivatives of fitted curves, smooth using Epanechnikov kernel
fittedCurcesP1 <- fitted(FPCAsparse, derOptns=list(p = 1, kernelType = 'epan'))

## ----eval=TRUE, echo=TRUE------------------------------------------------
  # load data
  data(medfly25)

  # Turn the original data into a list of paired amplitude and timing lists
  Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs)
  fpcaObjFlies <- FPCA(Flies$Ly, Flies$Lt, list(plot = TRUE, methodMuCovEst = 'smooth', userBwCov = 2))

## ----eval=TRUE, echo=TRUE------------------------------------------------
require('ks')
par(mfrow=c(1,2))
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), main = 'K = 11', pch = 4); grid()
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), K = 3, main = 'K = 3', pch = 4) ; grid()

## ----eval=TRUE, echo=TRUE------------------------------------------------
par(mfrow=c(1,2))
  CreateOutliersPlot(fpcaObjFlies, optns = list(K = 3, variant = 'KDE'))
  CreateFuncBoxPlot(fpcaObjFlies, xlab = 'Days', ylab = '# of eggs laid', optns = list(K =3, variant='bagplot'))

## ----eval=TRUE, echo=TRUE------------------------------------------------
par(mfrow=c(1,2))
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), K = 3, main = 'K = 3', showObs = FALSE) ; grid()
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), K = 3, main = 'K = 3', showObs = FALSE, derOptns = list(p = 1, bw = 1.01 , kernelType = 'epan') ) ; grid()

## ----eval=TRUE, echo=TRUE------------------------------------------------
fpcaObjFlies79 <- FPCA(Flies$Ly, Flies$Lt, list(nRegGrid = 79, methodMuCovEst = 'smooth', userBwCov = 2)) # Use 79 equidistant points for the support
CreateBWPlot(fpcaObjFlies79 , derOptns = list(p = 1, bw = 2.0 , kernelType = 'rect') )

## ----eval=TRUE, echo=TRUE------------------------------------------------
A <- FClust(Flies$Ly, Flies$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90), k = 2)
# The Neg-Entropy Criterion can be found as: A$clusterObj@bestResult@criterionValue 
CreatePathPlot( fpcaObjFlies, K=2, showObs=FALSE, lty=1, col= A$cluster, xlab = 'Days', ylab = '# of eggs laid')
grid()

