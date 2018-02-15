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
  cRates <- c( sum(trueClusters != C$cluster), sum(trueClusters != ifelse(D$cluster==1, 2, 1)))/N # 0.9677 & 0.9355 

  expect_gt( cRates[2], 0.935) # kCFC
  expect_gt( cRates[1], 0.967) # Rmixmod

  load(system.file('data', 'medfly25.RData', package='fdapace'))
  Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs) 
   for (i in 1:3) {
     set.seed(i)
     A <- FClust(Flies$Ly, Flies$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90), seed=i)
     # B <- FClust(Flies$Ly, Flies$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90), k = 2, seed=i)
   }
})

