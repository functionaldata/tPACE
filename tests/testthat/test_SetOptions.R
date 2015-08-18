 cat("\nTests for 'SetOptions'")

library(testthat)


 optns = list(bwmu = NULL, bwmuMethod = NULL, bwuserCov = NULL, bwuserCovGcv = NULL,
          ntest1 = NULL,  selectionMethod = NULL, FVEthreshold = NULL,
          maxK = NULL, dataType = NULL, error = NULL, nRegGrid = NULL, 
          methodXi = NULL, shrink = NULL, kernel = NULL, 
          numComponents = NULL, diagnosticsPlot = NULL,
          numBins = NULL, yname = NULL, rho = NULL, rotationCut = NULL,
          verbose = NULL, userMu = NULL, userCov = NULL, methodMu = NULL, methodCov = NULL, 
          outPercent = NULL, useBinnedData = NULL, rotationCut = NULL)
 
 test_that("SetOptions test: optns$method = NULL 'IN' for Dense case, 'CE' for Sparse case, 'CE' for other cases with warning", { 
  expect_equal(SetOptions(list(c(1,3,5), c(2,4)),list(c(1,3,5), c(2,4)), optns)$methodXi, 'CE') 
  expect_equal(SetOptions(list(c(1,2,3), c(1,2,3)),list(c(1,2,3), c(1,2,3)), optns)$methodXi, 'IN') 
  expect_equal(SetOptions(list(c(1,2,3,4,5), c(1,2,3,4)),list(c(1,2,3,4,5), c(1,2,3,4)), optns)$methodXi, 'IN') # DenseWithMV so should be 'IN' not 'CE'?
})
