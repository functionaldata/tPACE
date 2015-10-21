devtools::load_all()
library(testthat)

xcov <- (read.csv('data/xcovForGetEigens1.csv', header=FALSE))
xcov <-  do.call(rbind, xcov)
regGrid<- seq(0,10,length.out=51)

tmp <- no_FVE(xcov, 1, returnEVec=TRUE)
test_that('no_FVE interfaces correctly', {
  expect_equal( sum(tmp$lambda), 14.316031705480098 * 5)
 # expect_equal(length(tmp$lambda), length(tmp$FVE)) # Obs. check
  expect_equal(length(tmp$lambda), ncol(tmp$eVec))
})

