devtools::load_all(); library(testthat)

test_that('the growth example works.', {


         A <- read.table('data/growth.dat')
         B <- MakeFPCAInputs( IDs = A[,1], tVec = A$V3, yVec = A$V4) 
         C <- FClust(B$Ly, B$Lt, k = 2, cmethod = 'Rmixmod')
         D <- FClust(B$Ly, B$Lt, k = 2, cmethod = 'kCFC')
         trueClusters <-  A$V2[!duplicated(A$V1)]
         N = length(trueClusters)
         cRates <- c( sum(trueClusters != C$cluster), sum(trueClusters != ifelse(D$cluster==1, 2, 1)))/N # 0.9677 & 0.9355 

         expect_more_than( cRates[2], 0.935) # kCFC
         expect_more_than( cRates[1], 0.967) # Rmixmod
  
})

