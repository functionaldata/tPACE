 cat("\nTests for 'GetSmoothedMeanCurve.R'")

load('data/dataGeneratedByExampleSeed123.RData')

test_that("basic that the Epan. kernel gives the same results as MATLAB", {

  p = CreateOptions(kernel='epan')
  optns = SetOptions(y,t,p)
  out1 = sort(unique( c(unlist(t), optns$newdata)));
  out21 = seq(min(out1), max(out1),length.out = optns$ngrid);
  smcObj = GetSmoothedMeanCurve(y, t, out1, out21, p)
  expect_equal( sum(smcObj$mu) , 1.176558873333339e+02,tolerance = 1e-13, scale = 1 ) 

 } )

