cat("\nTests for 'gder'")

test_that("basic vector arguments do not return any errors / match the MATLAB output ", { 
  expect_equal( gder(c(1.0,2.0,3.0,4.0)), gder(c(1.0,2.0,3.0,4.0), family = 'poisson')) 
  expect_equal( gder(c(1.0,2.0,3.0,4.0), family = 'poisson'),  c(2.718281828459046, 7.389056098930650, 20.085536923187668, 54.598150033144236)) 
  expect_equal( gder(c(1.0,2.0,3.0,4.0), family = 'binomial'), c(0.196611933241482, 0.104993585403507,  0.045176659730912,  0.017662706213291)) 
  expect_equal( gder(c(0.1,0.2,0.3,0.4), family = 'binomial'), c(0.249376040192892, 0.247516572711860,  0.244458311690746,  0.240260745741529))
})
