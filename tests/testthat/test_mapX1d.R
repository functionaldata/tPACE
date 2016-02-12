cat("\ntests for 'MapX1D'")

test_that("basic arguments do not return any errors ", {
 xn = c(1:4,16)
 y = matrix(1:30, 15,2)
 x = c(1:14,16)
 expect_equal(  MapX1D(x,y,xn), matrix(c(1,2,3,4,15,16,17,18,19,30), 5,2))
 expect_equal(  MapX1D(1:14, seq(0,1,length.out=14),1:4 ),  seq(0,1,length.out=14)[1:4])                
}
)

