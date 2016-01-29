devtools::load_all()
library(testthat)
set.seed(1)
n <- 500
p <- 100
pts <- seq(0, 1, length.out=p)
# samp <- matrix(rnorm(n), n, p)
samp <- wiener(n, pts)
samp <- samp + rnorm(n * p, sd=0.01) + matrix(pts, n, p, byrow=TRUE)
# matplot(t(samp[1:10, ]), type='l')
samp <- sparsify(samp, pts, p)
res <- FPCA(samp$yList, samp$tList, 
            list(dataType='Dense', error=FALSE, lean=TRUE))
# fittedY <- fitted(res, derOptns=list(p=2))
# matplot(t(fittedY[1:10, ]), type='l', ylim=c(0, 2))


resDer <- deriv(res, list(p=1))
# plot(resDer$mu, type='l', ylim=c(-1, 1))
# plot(resDer$muDer, type='l', ylim=c(-2, 2))
# plot(resDer$phi[, 1], type='l', ylim=c(-2, 2))
# plot(resDer$phiDer[, 1], type='l', ylim=c(-2, 2))
test_that('first derivative is fine', {
  expect_equal(mean(resDer$muDer), 1, tolerance=0.1)
})
