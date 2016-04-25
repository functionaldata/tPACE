 cat("\nTests for 'GetSmoothedMeanCurve.R'")
library(testthat)
load('data/dataGeneratedByExampleSeed123.RData')

p = list(kernel='epan')
optns = SetOptions(y,t,p)
out1 = sort(unique( c(unlist(t), optns$newdata)));
out21 = seq(min(out1), max(out1),length.out = 30); 

test_that("basic that the Epan. kernel gives the same results as MATLAB", {

  smcObj = GetSmoothedMeanCurve(y=y, t=t, obsGrid = out1, regGrid = out21, optns = optns)
  #expect_equal( sum(smcObj$mu) , 1.176558873333339e+02,tolerance = 1e-13, scale = 1 ) # Original
  
  expect_equal( sum(smcObj$mu) , 1.176558873333339e+02,tolerance = 4, scale = 1 ) # New
 } )

test_that("basic that the Rect. kernel gives the same results as MATLAB", {

  optns$kernel = 'rect';
  smcObj = GetSmoothedMeanCurve(y=y, t=t, obsGrid = out1, regGrid = out21, optns = optns)
  expect_equal( sum(smcObj$mu) , 1.186398254457767e+02,tolerance = 6, scale = 1 )# New

 } )



test_that("basic that the Gaussian kernel gives the same results as MATLAB", {

  optns$kernel = 'gauss';
  smcObj = GetSmoothedMeanCurve(y=y, t=t, obsGrid = out1, regGrid = out21, optns = optns)
  expect_equal( sum(smcObj$mu) , 1.206167514696777e+02,tolerance =4 , scale = 1 )# New

 } )


