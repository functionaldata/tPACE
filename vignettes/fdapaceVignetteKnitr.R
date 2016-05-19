## ----include=FALSE-------------------------------------------------------
library("knitr") 
opts_chunk$set(warning=FALSE, fig.width=5, fig.height=5, fig.align='center', cache=TRUE, background='white', highlight=FALSE)
# Set hooks
render_listings()

## ----eval=FALSE----------------------------------------------------------
#  FPCAobj <- FPCA(Ly=yList, Lt=tList)

## ----message=FALSE-------------------------------------------------------
  #library(fdapace)
  devtools::load_all('.') # So we use the new interfaces

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

## ----denseSim, fig.show='hide'-------------------------------------------
  L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M),tVec=rep(s,N), t(yTrue))
  FPCAdense <- FPCA(L3$Ly, L3$Lt)   
  
  # Make a basic diagnostics plot
  plot(FPCAdense)
  
  # Find the standard deviation associated with each component
  sqrt(FPCAdense$lambda)

## ----sparseSim, fig.show='hide'------------------------------------------
  # Create sparse sample  
  # Each subject has one to five readings (median: 3)
  set.seed(123)
  ySparse <- Sparsify(yTrue, s, sparsity = c(1:5))
    
  # Give your sample a bit of noise 
  ySparse$yNoisy <- lapply(ySparse$Ly, function(x) x + 0.5*rnorm(length(x))) 
  
  
  # Do FPCA on this sparse sample
  # Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
  # Smoothing is the main computational cost behind sparse FPCA
  FPCAsparse <- FPCA(ySparse$yNoisy, ySparse$Lt, list(diagnosticsPlot = TRUE))  
  

## ------------------------------------------------------------------------
 FPCAsparseMuBW5 <- FPCA(ySparse$yNoisy, ySparse$Lt, optns= list(userBwMu = 5)) 

## ----bwChoice, fig.show='hold', out.width='.45\\linewidth'---------------
CreatePathPlot( FPCAsparse, subset = 1:3, main = "GCV bandwidth", pch = 16)
CreatePathPlot( FPCAsparseMuBW5, subset = 1:3, main = "User-defined bandwidth", pch = 16)

## ------------------------------------------------------------------------
 FPCAsparseRect <- FPCA(ySparse$yNoisy, ySparse$Lt, optns = list(kernel = 'rect')) # Use rectangular kernel

## ----eval=FALSE----------------------------------------------------------
#  SelectK( FPCAsparse, criterion = 'FVE', FVEthreshold = 0.95) # k = 2
#  SelectK( FPCAsparse, criterion = 'AIC') # k = 2

## ------------------------------------------------------------------------
fittedCurvesP0 <- fitted(FPCAsparse) # equivalent: fitted(FPCAsparse, derOptns=list(p = 0));
# Get first order derivatives of fitted curves, smooth using Epanechnikov kernel
fittedCurcesP1 <- fitted(FPCAsparse, derOptns=list(p = 1, kernelType = 'epan')) 

## ----medfly, fig.show='hide'---------------------------------------------
  # load data
  data(medfly25)

  # Turn the original data into a list of paired amplitude and timing lists
  Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs) 
  fpcaObjFlies <- FPCA(Flies$Ly, Flies$Lt, list(diagnosticsPlot = TRUE, methodMuCovEst = 'smooth', userBwCov = 2)) 

## ----k, fig.show='hold', out.width='.45\\linewidth'----------------------
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), main = 'k = 11', pch = 4); grid()
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), k = 3, main = 'k = 3', pch = 4) ; grid()  

## ----outlier, fig.show='hold', out.width='.4\\linewidth'-----------------
  CreateOutliersPlot(fpcaObjFlies, optns = list(k = 3, variant = 'KDE'))
  CreateFuncBoxPlot(fpcaObjFlies, xlab = 'Days', ylab = '# of eggs laid', optns = list(k =3, variant='bagplot'))

## ----deriv, fig.show='hold', out.width='.45\\linewidth'------------------
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), k = 3, main = 'k = 3', showObs = FALSE) ; grid() 
  CreatePathPlot(fpcaObjFlies, subset = c(3,5,135), k = 3, main = 'k = 3', showObs = FALSE, derOptns = list(p = 1, bw = 1.01 , kernelType = 'epan') ) ; grid() 

## ----bwPlot, fig.width=7, fig.height=3-----------------------------------
fpcaObjFlies79 <- FPCA(Flies$Ly, Flies$Lt, list(nRegGrid = 79, methodMuCovEst = 'smooth', userBwCov = 2)) # Use 79 equidistant points for the support
CreateBWPlot(fpcaObjFlies79 , derOptns = list(p = 1, bw = 2.0 , kernelType = 'rect') )

