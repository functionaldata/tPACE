# devtools::load_all()
library(testthat)

test_that('MakeFPCAInputs works', {
  set.seed(1)
  NN = 32
  LyOrig = sapply( 1:NN, function(x)  rnorm(x), simplify=FALSE)
  LtOrig = sapply( 1:NN, function(x)  sample(1:100,x), simplify=FALSE)

  LyOrigVec = as.vector(unlist(LyOrig))
  LtOrigVec = as.vector(unlist(LtOrig))

  LusersOrig = sapply( 1:NN, function(x) rep( paste0(collapse='', 'user', x), x) )
  LusersOrigVec = as.vector(unlist(LusersOrig))
  ZZ = MakeFPCAInputs(IDs= LusersOrigVec, tVec=LtOrigVec,yVec= LyOrigVec)
  # BB = FPCA(Ly= ZZ$Ly, Lt= ZZ$Lt) # This errs!
  ZZs = MakeFPCAInputs(IDs= LusersOrigVec, tVec=LtOrigVec,yVec= LyOrigVec, beSorted=TRUE)
  # CC = FPCA(Ly= ZZs$Ly, Lt= ZZs$Lt)

  expect_s3_class(FPCA(Ly= ZZs$Ly, Lt= ZZs$Lt), 'FPCA')
  expect_error( FPCA(Ly = ZZ$Ly, Lt = ZZ$Lt), "Each vector in t should be in ascending order" )

})
