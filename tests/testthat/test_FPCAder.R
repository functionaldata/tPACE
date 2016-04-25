devtools::load_all()
library(testthat)

test_that('FPCAder correct derivatives of mean for dense case', {
  set.seed(1)
  n <- 500
  p <- 300
  pts <- seq(0, 1, length.out=p)
  samp <- Wiener(n, pts)
  samp <- samp + rnorm(n * p, sd=0.01) + matrix(pts, n, p, byrow=TRUE)
  spSamp <- Sparsify(samp, pts, p) # This just puts the sample in list.
  fpcaObj <- FPCA(spSamp$Ly, spSamp$Lt, list(dataType='Dense' ))
  fpcaObjDer <- FPCAder(fpcaObj, list(p=1))
    
  expect_equal(median(fpcaObjDer$muDer), 1, tolerance=0.1, scale = 1)
})




test_that('FPCAder correct derivatives of mean for sparse case',{
  set.seed(1)
  n <- 333 
  pts <- seq(0, 1, by=0.01)
  mu <- 0:(length(pts)-1) / 50
  phi =  CreateBasis(K=4, type='fourier',  pts= pts)
  samp1 <- t(mu + phi %*% matrix( rnorm(n*4, mean=0, sd=c(1,.4, 0.01, 0.001)), nrow=4))
  samp2 <- Sparsify(samp1, pts, 8)
  fpcaObj = FPCA(samp2$Lt, Ly= samp2$Ly)
  fpcaObjDer <- FPCAder(fpcaObj, list(p=1))
  expect_equal(median(fpcaObjDer$muDer), 2, tolerance=0.1, scale = 1) 
  })
