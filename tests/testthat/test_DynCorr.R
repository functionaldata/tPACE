#devtools::load_all()
library(testthat)

#positive test 1
test_that('when y_i[t] - mean(y[t]) = 2(x_i[t]-mean(x[t])), expect correlation 1', {
  library(MASS)
  set.seed(10)
  n=200 
  t=seq(0,1,length.out=100)
  
  mu_quad_x=8*t^2-4*t+5
  mu_quad_y=8*t^2-12*t+6
  
  fun=rbind(rep(1,length(t)),-t,t^2)
  z1=mvrnorm(n,rep(0,3),diag(c(2,16/3,4)))   # covariance matrix of random effects
  x1_quad=y1_quad=x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))
  
  for (i in 1:n){
    x1_quad[i,]=mu_quad_x+z1[i,]%*%fun
    y1_quad[i,]=mu_quad_y+2*z1[i,]%*%fun
  }
  
  dyn1_quad=DynCorr(x1_quad,y1_quad,t)
  expect_equal(min(dyn1_quad),1)
})

#positive test 2
test_that('first generate y_i[t] - mean(y[t]) = 2(x_i[t]-mean(x[t])) then add noise to y and x, expect correlation close to 1', {
  library(MASS)
  set.seed(10)

  n=200 
  t=seq(0,1,length.out=100)

  mu_quad_x=8*t^2-4*t+5
  mu_quad_y=8*t^2-12*t+6

  fun=rbind(rep(1,length(t)),-t,t^2)
  z1=mvrnorm(n,rep(0,3),diag(c(2,16/3,4)))   # covariance matrix of random effects
  x1_quad=y1_quad=x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))

  for (i in 1:n){
    x1_quad[i,]=mu_quad_x+z1[i,]%*%fun+rnorm(length(t),0,0.01)
    y1_quad[i,]=mu_quad_y+2*z1[i,]%*%fun+rnorm(length(t),0,0.01)
  }

  dyn1_quad=DynCorr(x1_quad,y1_quad,t)
  expect_equal(max(1-dyn1_quad) < 0.1,TRUE)
})

#positive test 3
test_that('setting same as test #1 but with 1/5 missings in each observation, expect result close to 1', {
  library(MASS)
  set.seed(10)
  
  n=200 
  t=seq(0,1,length.out=100)
  
  mu_quad_x=8*t^2-4*t+5
  mu_quad_y=8*t^2-12*t+6
  
  fun=rbind(rep(1,length(t)),-t,t^2)
  z1=mvrnorm(n,rep(0,3),diag(c(2,16/3,4)))   # covariance matrix of random effects
  x1_quad=y1_quad=x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))
  
  for (i in 1:n){
    x1_quad[i,]=mu_quad_x+z1[i,]%*%fun
    y1_quad[i,]=mu_quad_y+2*z1[i,]%*%fun
  }
  
  for(i in 1:n){
    ms_x = sample(t,20)
    ms_y = sample(t,20)
    x1_quad[i,ms_x] = NA
    y1_quad[i,ms_y] = NA
  }
  
  dyn1_quad=DynCorr(x1_quad,y1_quad,t)
  expect_equal(max(1-dyn1_quad) < 0.1,TRUE)
})

#negative test
expect_error({
  set.seed(10)
  n=200 
  t=seq(0,1,length.out=100)
  
  x = Wiener(n,t)
  y= Wiener(n+1,t)
  
  dyn=DynCorr(x,y,t)
})