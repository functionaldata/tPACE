#devtools::load_all()
library(testthat)

#positive test 1
test_that('test estimates when x,y are independent', {
  set.seed(10)
  n = 20000
  m = 50
  ei = rnorm(n)
  y=1:n
  x=list()
  t_x=list()
  for(i in 1:n){
    t_x = c(t_x,list(0:m/m))
    x = c(x,list(ei[i]*array(1,m+1)))
    y[i] = rnorm(1)
  }
  
  outQ = c(0.1,0.25,0.5,0.75,0.9,0.95)
  isNewsub = rep(0,n)
  qtreg = FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewsub)
  
  expect_equal(max(abs(qtreg$pred_quantile[,1]-qnorm(0.1)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[,2]-qnorm(0.25)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[,3]-qnorm(0.5)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[,4]-qnorm(0.75)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[,5]-qnorm(0.9)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[,6]-qnorm(0.95)))<0.1,T)
})

#positive test 2
test_that('test estimates when x only takes 2 values', {
  set.seed(10)
  n = 2000
  m = 50
  ei = c(rep(0,n/2),rep(1,n/2))
  y=1:n
  x=list()
  t_x=list()
  for(i in 1:n){
    t_x = c(t_x,list(0:m/m))
    x = c(x,list(ei[i]*array(1,m+1)))
    if(ei[i]==0){
      y[i] = rnorm(1)
    }else{
      y[i] = runif(1,0,1)
    }
  }
  
  outQ = c(0.1,0.25,0.5,0.75,0.9)
  isNewsub = rep(0,n)
  qtreg = FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewsub)
  
  expect_equal(max(abs(qtreg$pred_quantile[ei==0,1]-qnorm(0.1)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==0,2]-qnorm(0.25)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==0,3]-qnorm(0.5)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==0,4]-qnorm(0.75)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==0,5]-qnorm(0.9)))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==1,1]-0.1))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==1,2]-0.25))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==1,3]-0.5))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==1,4]-0.75))<0.1,T)
  expect_equal(max(abs(qtreg$pred_quantile[ei==1,5]-0.9))<0.1,T)
})


# negative test 1
test_that('test when output quantiles are not in [0,1]', {
  n = 200
  npred = 50
  m = 50
  xi <- Wiener(n, 0:m/m)
  y=1:n
  x=list()
  t_x=list()
  for(i in 1:n){
    t_x = c(t_x,list(0:m/m))
    x = c(x,list(xi[i,]))
    y[i] = 5*rnorm(1)+2*sum(xi[i,])
  }
  
  outQ = c(-0.1,0.25,0.5,0.75,0.9,0.95)
  isNewsub = c(rep(0,n-npred),rep(1,npred))
  expect_error(FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewsub) )
  #qtreg = FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewsub) 
})

# negative test 2
test_that('test when y and Lx does not have same length', {
  n = 200
  npred = 50
  m = 50
  xi <- Wiener(n, 0:m/m)
  y=1:(n+1)
  x=list()
  t_x=list()
  for(i in 1:n){
    t_x = c(t_x,list(0:m/m))
    x = c(x,list(xi[i,]))
    y[i] = 5*rnorm(1)+2*sum(xi[i,])
  }
  
  outQ = c(0.1,0.25,0.5,0.75,0.9,0.95)
  isNewsub = c(rep(0,n-npred),rep(1,npred))
  expect_error(FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewsub) )
  #qtreg = FPCquantile(x, t_x, y, outQ,optns_x = NULL,isNewsub) 
})
