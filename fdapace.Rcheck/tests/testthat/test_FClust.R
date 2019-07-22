# devtools::load_all()
library(testthat)

test_that('the growth example works.', {


  A <- read.table(system.file('testdata', 'growth.dat',
                             package='fdapace'))
  B <- MakeFPCAInputs( IDs = A[,1], tVec = A$V3, yVec = A$V4) 
  C <- FClust(B$Ly, B$Lt, k = 2, cmethod = 'EMCluster')
  D <- FClust(B$Ly, B$Lt, k = 2, cmethod = 'kCFC')
  trueClusters <-  A$V2[!duplicated(A$V1)]
  N = length(trueClusters)
  cRates <- c( sum(trueClusters == C$cluster), sum(trueClusters == C$cluster) )/N # 0.9677 & 0.9355 
  cRates <- sapply(cRates, function(x) ifelse(x < 0.5, 1- x, x))

  expect_gt( cRates[2], 0.935) # kCFC
  expect_gt( cRates[1], 0.967) # Rmixmod

  load(system.file('data', 'medfly25.RData', package='fdapace'))
  Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs) 
   for (i in 1:3) {
     set.seed(i)
     A <- FClust(Flies$Ly, Flies$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90))
     # B <- FClust(Flies$Ly, Flies$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90), k = 2, seed=i)
   }
})

test_that('the k-means initialisation error occurs normally.', {
  
  set.seed(1)
  n <- 100
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  sampTrue <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE)
  samp <- sampTrue + rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- MakeFPCAInputs(tVec=pts, yVec=samp)
  expect_error(FClust(tmp$Ly, tmp$Lt, k = 14, cmethod = "kCFC",  optnsFPCA = list(userBwCov= 2, FVEthreshold = 0.90)) )
  
})
