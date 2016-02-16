
myEps <- .Machine$double.eps;

test_that(" basic 0/1 combination matches MATLAB output ", { 
  expect_equal( CreateTrueMean(1,0), 0)
  expect_equal( CreateTrueMean(0,1), 0)
  expect_equal( CreateTrueMean(0,0), 0)
  expect_equal( CreateTrueMean(1,1), 1.841470984807897, tolerance = 2*myEps, scale = 1)
})

test_that(" 'multiple t, single p'-case matches MATLAB output ", { 
  expect_equal( CreateTrueMean( c(0.1, 1.3), 3), c(0.199833416646828, 2.263558185417193), tolerance = 2*myEps, scale = 1 ) 
  expect_equal( CreateTrueMean( c(0.1, 1.3), 1), c(0.199833416646828, 0) , tolerance = 2*myEps, scale = 1 ) 
})

test_that(" 'single t, multiple p'-case matches MATLAB output ", { 
  expect_equal( CreateTrueMean( 1, c(0.1, 1.3)), 0 ) 
  expect_equal( CreateTrueMean( 3, c(0.1, 1.3)), c(0, 0)  ) 
})

test_that(" 'multiple t, multiple p'-case matches MATLAB output ", { 
  expect_equal( CreateTrueMean( c(3, 1), c(0.1, 1.3)), c(0, 1.841470984807897), tolerance = 2*myEps, scale = 1)  
  expect_equal( CreateTrueMean( c(1, 3), c(0.1, 1.3)), c(0, 0)  , tolerance = 2*myEps, scale = 1) 
  expect_equal( CreateTrueMean( c(0.1, 1.3), c(3, 1)), c(0.199833416646828, 0), tolerance = 2*myEps, scale = 1)  
  expect_equal( CreateTrueMean( c(0.1, 1.3), c(1, 3)), c(0.199833416646828, 2.263558185417193), tolerance = 2*myEps, scale = 1)
  expect_equal( CreateTrueMean( c(1.1, 0.3), c(3, 1)), c(1.991207360061436, 0.595520206661340), tolerance = 2*myEps, scale = 1)  
  expect_equal( CreateTrueMean( c(1.1, 0.3), c(1, 3)), c(0, 0.595520206661340), tolerance = 2*myEps, scale = 1)  
})

# cat("Done")