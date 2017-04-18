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
  ZZs = MakeFPCAInputs(IDs= LusersOrigVec, tVec=LtOrigVec,yVec= LyOrigVec, sort=TRUE)
  # CC = FPCA(Ly= ZZs$Ly, Lt= ZZs$Lt)
  
  expect_s3_class(FPCA(Ly= ZZs$Ly, Lt= ZZs$Lt), 'FPCA')
  expect_error( FPCA(Ly = ZZ$Ly, Lt = ZZ$Lt), "Each vector in t should be in ascending order" )
  
})

test_that("basic arguments give reasonable output ", { 
  
  IDs = c('a','a','b','b','U34')
  tVec = 1:5;
  yVec = cos(tVec);
  B = MakeFPCAInputs(IDs= IDs , tVec=tVec, yVec=yVec)
  
  expect_equal( unlist(B$Lt), tVec, tolerance = 2*.Machine$double.eps, scale = 1) 
  expect_equal( B$Ly[[2]], cos(c(3,4)), tolerance = 2*.Machine$double.eps, scale = 1) 
  expect_true( (length(B$Lid) == length(B$Ly)) && (length(B$Ly) == length(B$Lt)) )
  expect_true( B$Lid[[3]] == IDs[5]  )
})

test_that("basic arguments give reasonable output when number of measurement points is equal ", {
  
  IDs = rep(1:3,each=3); 
  tVec = rep(c(0,2,5),3); 
  yVec = 10:19;
  
  B = MakeFPCAInputs(IDs= IDs , tVec=tVec, yVec=yVec)
  
  expect_equal( unlist(B$Lt), tVec, tolerance = 2*.Machine$double.eps, scale = 1)
  expect_equal( B$Ly[[2]], c(13,14,15), tolerance = 2*.Machine$double.eps, scale = 1)
  expect_true( (length(B$Lid) == length(B$Ly)) && (length(B$Ly) == length(B$Lt)) )
  expect_true( B$Lid[[3]] == IDs[9]  )
})


