cat("\nTests for 'ginver'")

test_that("basic vector arguments do not return any errors / match the MATLAB output ", { 
  expect_equal( ginver(c(1.0,2.0,3.0,4.0)), ginver(c(1.0,2.0,3.0,4.0), family = 'poisson')) 
  expect_equal( ginver(c(1.0,2.0,3.0,4.0), family = 'poisson'),  c( 0, 0.693147180559945, 1.098612288668110, 1.386294361119891)) 
  expect_equal( ginver(c(1.0,2.0,3.0,4.0), family = 'binomial'), rep( 4.993239250550510, times= 4) ) 
  expect_equal( ginver(c(0.1,0.2,0.3,0.4), family = 'binomial'), c(-2.197224577336219, -1.386294361119891, -0.847297860387204, -0.405465108108164))
})

# cat("Done")
