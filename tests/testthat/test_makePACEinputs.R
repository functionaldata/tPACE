 cat("\nTests for 'makePACEinputs'")


test_that("basic arguments give reasonable output ", { 

  IDs = c('a','a','b','b','U34')
  tVec = 1:5;
  yVec = cos(tvec);
  B = makePACEinputs(IDs= IDs , tVec=tvec, yVec=yvec)

  expect_equal( unlist(B$Lt), tvec, tolerance = 2*.Machine$double.eps, scale = 1) 
  expect_equal( B$Ly[[2]], cos(c(3,4)), tolerance = 2*.Machine$double.eps, scale = 1) 
  expect_true( (length(B$Lid) == length(B$Ly)) && (length(B$Ly) == length(B$Lt)) )
  expect_true( B$Lid[[3]] == IDs[5]  )
})
# cat("Done")
