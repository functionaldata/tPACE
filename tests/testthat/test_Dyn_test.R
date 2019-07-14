#devtools::load_all()
library(testthat)

#positive test 1
test_that('one sample test, when y_i[t] - mean(y[t]) = 2(x_i[t]-mean(x[t])), expect p-value very small', {
  library(MASS)
  set.seed(10)
  n=200 
  t=seq(0,1,length.out=100)
  
  mu_quad_x=8*t^2-4*t+5
  mu_quad_y=8*t^2-12*t+6
  
  fun=rbind(rep(1,length(t)),-t,t^2)
  z1=mvrnorm(n,rep(0,3),diag(c(2,16/3,4)))   # covariance matrix of random effects
  x1_quad=y1_quad=matrix(0,nrow=n,ncol=length(t))
  
  for (i in 1:n){
    x1_quad[i,]=mu_quad_x+z1[i,]%*%fun
    y1_quad[i,]=mu_quad_y+2*z1[i,]%*%fun
  }
  
  bt_dc=Dyn_test(x1_quad,y1_quad,t)
  expect_equal(bt_dc$pval < 0.01,TRUE)
})

#positive test 2
test_that('one sample test, setting same as test #1 but with 1/5 missings in each observation, expect p-value very small', {
  library(MASS)
  set.seed(10)
  
  n=200 
  t=seq(0,1,length.out=100)
  
  mu_quad_x=8*t^2-4*t+5
  mu_quad_y=8*t^2-12*t+6
  
  fun=rbind(rep(1,length(t)),-t,t^2)
  z1=mvrnorm(n,rep(0,3),diag(c(2,16/3,4)))   # covariance matrix of random effects
  x1_quad=y1_quad=matrix(0,nrow=n,ncol=length(t))
  
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
  
  bt_dc=Dyn_test(x1_quad,y1_quad,t)
  expect_equal(bt_dc$pval < 0.01,TRUE)
})

#positive test 3
test_that('two sample test, the first sample has the same setting as test #2, the second sample is the first sample with noise, expect p-value larger than 0.2', {
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
    x1_quad_error[i,]=x1_quad[i,]+rnorm(length(t),0,0.1)
    y1_quad_error[i,]=y1_quad[i,]+rnorm(length(t),0,0.1) 
  }
  
  for(i in 1:n){
    ms_x = sample(t,20)
    ms_y = sample(t,20)
    x1_quad[i,ms_x] = NA
    y1_quad[i,ms_y] = NA
    x1_quad_error[i,ms_x] = NA
    y1_quad_error[i,ms_y] = NA
  }
  
  bt_dc=Dyn_test(x1_quad,y1_quad,t,x1_quad_error,y1_quad_error,t)
  expect_equal(bt_dc$pval > 0.2,TRUE)
})
