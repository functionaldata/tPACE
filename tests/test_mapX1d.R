cat("\ntests for 'mapX1d'")

test_that("basic arguments do not return any errors ", {
 xn = c(1:4,16)
 y = matrix(1:30, 15,2)
 x = c(1:14,16)
 expect_equal(  mapX1d(x,y,xn), matrix(c(1,2,3,4,15,16,17,18,19,30), 5,2))
 expect_equal(  mapX1d(1:14, seq(0,1,14),1:4 ),  seq(0,1,14)[1:4])
                
}
)

