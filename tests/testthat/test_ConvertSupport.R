devtools::load_all()

fromGrid <- seq(0, pi / 2, 0.1)
toGrid <- fromGrid + 0.001
toGrid[length(toGrid)] <- toGrid[length(toGrid)] - 0.002
mu <- sin(fromGrid)
phi <- cbind(sin(fromGrid), cos(fromGrid))
phi1 <- matrix(phi[, 1], ncol=1)
Cov <- tcrossprod(phi)

test_that('ConvertSupport works', {
  expect_equal(mu, ConvertSupport(fromGrid, toGrid, mu=mu), tolerance=2e-3)
  expect_equal(phi, ConvertSupport(fromGrid, toGrid, phi=phi), tolerance=2e-3)
  expect_equal(Cov, ConvertSupport(fromGrid, toGrid, Cov=Cov), tolerance=1e-3)
})
