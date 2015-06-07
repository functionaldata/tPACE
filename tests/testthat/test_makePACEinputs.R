 cat("\nTests for 'makePACEinputs'")


test_that("basic arguments give reasonable output ", { 

  IDs = c('a','a','b','b','U34')
  tVec = 1:5;
  yVec = cos(tVec);
  B = makePACEinputs(IDs= IDs , tVec=tVec, yVec=yVec)

  expect_equal( unlist(B$Lt), tVec, tolerance = 2*.Machine$double.eps, scale = 1) 
  expect_equal( B$Ly[[2]], cos(c(3,4)), tolerance = 2*.Machine$double.eps, scale = 1) 
  expect_true( (length(B$Lid) == length(B$Ly)) && (length(B$Ly) == length(B$Lt)) )
  expect_true( B$Lid[[3]] == IDs[5]  )
})

test_that("basic arguments give reasonable output when number of measurement points is equal ", {

  IDs = rep(1:3,each=3); 
  tVec = rep(c(0,2,5),3); 
  yVec = 10:19;

  B = makePACEinputs(IDs= IDs , tVec=tVec, yVec=yVec)

  expect_equal( unlist(B$Lt), tVec, tolerance = 2*.Machine$double.eps, scale = 1)
  expect_equal( B$Ly[[2]], c(13,14,15), tolerance = 2*.Machine$double.eps, scale = 1)
  expect_true( (length(B$Lid) == length(B$Ly)) && (length(B$Ly) == length(B$Lt)) )
  expect_true( B$Lid[[3]] == IDs[9]  )
})




# cat("Done")
