cat("\ntests for 'getEigens'")

myEps = .Machine$double.eps;
test_that("basic arguments match MATLAB output ", { 

 xcov <- (read.csv('../../data/xcovForGetEigens1.csv', header=FALSE))
 xcov <-  do.call(rbind, xcov)
 obsGrid<- seq(0,10,length.out=21)
 regGrid<- seq(0,10,length.out=51)
 noeig <- 23;
 AA <- getEigens(xcov, obsGrid=obsGrid, regGrid= regGrid, noeig=23)

  expect_equal( sum(AA$lambda), 14.316031705480098, tolerance = 100*myEps, scale = 1) 
  expect_equal( sum(AA$phi), -3.472755554844390, tolerance = 1e-10, scale = 1) 
  expect_equal( sum(AA$eigen),  -7.934784142851369 , tolerance = 1e-11, scale = 1)   
})
