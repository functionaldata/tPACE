library(testthat)
# devtools::load_all('.')

test_that('FuncCorrCent works', {
  set.seed(4)

  n <- 200
  nGridIn <- 50
  sparsity <- 1:5 # must have length > 1
  bw <- 0.1
  kern <- 'epan'
  T <- matrix(seq(0.5, 1, length.out=nGridIn))

  ## Corr(X(t), Y(t)) = 1/2
  A <- Wiener(n, T)
  B <- Wiener(n, T) 
  C <- Wiener(n, T) + matrix((1:nGridIn) , n, nGridIn, byrow=TRUE)
  X <- A + B
  Y <- A + C
  indEach <- lapply(1:n, function(x) sort(sample(nGridIn, sample(sparsity, 1))))
  tAll <- lapply(1:n, function(i) T[indEach[[i]]])
  Xsp <- lapply(1:n, function(i) X[i, indEach[[i]]])
  Ysp <- lapply(1:n, function(i) Y[i, indEach[[i]]])

  expect_equal(sapply(Xsp, length), sapply(Ysp, length))
  expect_equal(sapply(Xsp, length), sapply(tAll, length))

  # Perfect correlation case
  expect_equal(mean(FuncCorrCent(Xsp, Xsp, tAll, bw)[['corr']], na.rm=TRUE), 1)

  # Consistency
  expect_true(mean((FuncCorrCent(Xsp, Ysp, tAll, bw, kern)[['corr']] - 0.5)^2, na.rm=TRUE) < 1e-2)

  # Gauss and epan kernels are similar
  expect_equal(FuncCorrCent(Xsp, Ysp, tAll, bw, 'epan')[['corr']], FuncCorrCent(Xsp, Ysp, tAll, bw, 'gauss')[['corr']], 0.1)
})


