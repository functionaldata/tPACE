library(testthat)

load(system.file('testdata', 'dataGeneratedByExampleSeed123.RData', package='fdapace'))

test_that("sparse case", {
  expect_warning(res <- GetMeanCI(Lt=t, Ly=y, optns = list(nRegGrid = 30)))
  expect(all(res$CI$upper-res$CI$lower > 0), "upper bounds are not all greater than lower bounds.")
})

test_that("dense case", {
  n <- 30
  tgrid <- seq(0,1,length.out=21)
  phi1 <- function(t) sqrt(2)*sin(2*pi*t)
  phi2 <- function(t) sqrt(2)*sin(4*pi*t)
  Lt <- rep(list(tgrid), n)
  Ly <- lapply(1:n, function(i){
    tgrid + rnorm(1,0,2) * phi1(tgrid) + rnorm(1,0,0.5) * phi2(tgrid) + rnorm(1,0,0.01)
  })
  res <- GetMeanCI(Lt = Lt, Ly = Ly, level = 0.9)
  expect(all(res$CI$upper-res$CI$lower > 0),
         "upper bounds are not all greater than lower bounds.")
  res2 <- GetMeanCI(Lt = Lt, Ly = Ly, level = 0.95)
  expect(all(res2$CI$upper-res$CI$upper > 0) & all(res2$CI$lower - res$CI$lower < 0),
         "95% CIs are not wider than 90%.")
})
