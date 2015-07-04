
myEps <- .Machine$double.eps;
load('data/dataForGetRawCov.RData')
AA = GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',TRUE)  #Matches ML output
BB = GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output

test_that(" basic argument match MATLAB output ", { 
  expect_equal( sum(AA$indx) , 184, tolerance = 2*myEps, scale = 1)
  expect_equal( sum(AA$cxxn) , -7.416002855888680, tolerance = 1e-13, scale = 1)
  expect_equal( sum(AA$cyy) , 16.327874649330514, tolerance = 1e-13, scale = 1)
  expect_equal( sum(AA$tPairs) , 4.053285461728229e+02, tolerance = 1e-12, scale = 1)  
})

test_that(" basic argument match MATLAB output ", { 
  expect_equal( sum(BB$indx) , 298, tolerance = 2*myEps, scale = 1)
  expect_equal( sum(BB$cxxn) , 16.327874649330514, tolerance = 1e-13, scale = 1)
  expect_equal( sum(BB$cyy) , 16.327874649330514, tolerance = 1e-13, scale = 1)
  expect_equal( sum(BB$tPairs) , 6.330209554605514e+02, tolerance = 1e-12, scale = 1)  
})

y2 = list(1:10, 2:11)
t2 = list( 1:10, 1:10)

CC = GetRawCov(y2,t2, sort(unique(unlist(t2))), seq(1.5,10.5, length.out=10) ,'Dense',TRUE) #Matches ML output
DD = GetRawCov(y2,t2, sort(unique(unlist(t2))), seq(1.5,10.5, length.out=10) ,'Dense',FALSE) #Matches ML output
# DD = getRawCov(y2,t2, sort(unique(cell2mat(t2))), linspace(1.5,10.5, 10), 2, 0)


test_that(" basic argument match MATLAB output ", { 
  expect_equal( sum(CC$indx) , 0, tolerance = 2*myEps, scale = 1)
  expect_equal( sum(CC$cxxn) , 22.5, tolerance = 1e-13, scale = 1)
  expect_equal( sum(CC$cyy) , 25, tolerance = 1e-13, scale = 1)
  expect_equal( sum(CC$tPairs) , 990, tolerance = 1e-12, scale = 1)  
})


test_that(" basic argument match MATLAB output ", { 
  expect_equal( sum(DD$indx) , 0, tolerance = 2*myEps, scale = 1)
  expect_equal( sum(DD$cxxn) , 25, tolerance = 1e-13, scale = 1)
  expect_equal( sum(DD$cyy) , 25, tolerance = 1e-13, scale = 1)
  expect_equal( sum(DD$tPairs) , 1100, tolerance = 1e-12, scale = 1)  
})


