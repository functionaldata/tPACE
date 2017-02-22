# test the performance of Stringing
library(testthat)

# implement the first kind of simulation in stringing paper
# Fourier basis, specify signal-to-noise ratio (SNR), number of components used
# returns true underlying curves XTraj, and randomly ordered design matrix X
stringing_sim1 <- function(K = 4, SNR = 10, n = 50, p = 50){
  grid = seq(5/p,5,5/p)
  lambda = 8*(0.5^(1:K)) # exponential decay
  phi = cbind(-sqrt(0.2)*cos(0.2*pi*grid), sqrt(0.2)*sin(0.2*pi*grid),
              -sqrt(0.2)*cos(0.4*pi*grid), sqrt(0.2)*sin(0.4*pi*grid))
  scores = t(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(lambda)))
  Xtrue = t(phi %*% scores) # ith row corresponds to ith subject
  if(SNR == Inf){
    Xt = Xtrue
  } else {
    Xt = Xtrue + rnorm(n*p, mean=0, sd=mean(abs(Xtrue))/SNR)
  }
  rdmorder = sample(x = 1:p, size = p, replace = FALSE) # random columns indices for the original order
  X = Xt[,rdmorder]
  return(list(XTraj = Xt, X = X, order = rdmorder, n = n, p = p, RegGrid = grid))
}

test_that("Stringing works perfectly for simulated example without error using correlation metric",{
  set.seed(1)
  simdata = stringing_sim1(SNR = Inf)
  stringingfit = Stringing(simdata$X, disOptns = "correlation")
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  #plot(1:simdata$p, simdata$order[stringingfit$stringedPos], pch=18, 
  #     xlab="True Order", ylab="Stringed Order")
  #lines(1:simdata$p, simdata$order[stringingfit$stringedPos])
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_equal(1:simdata$p, simdata$order[stringingfit$stringedPos])
  expect_equal(ROE, 0)
})

test_that("Stringing works for simulated example with small error using correlation metric",{
  set.seed(1)
  simdata = stringing_sim1(SNR = 10)
  stringingfit = Stringing(simdata$X, disOptns = "correlation")
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_equal(1:simdata$p, simdata$order[stringingfit$stringedPos])
  expect_equal(ROE, 0)
})

test_that("Stringing works for simulated example with moderate error using correlation metric",{
  set.seed(2)
  simdata = stringing_sim1(SNR = 4)
  stringingfit = Stringing(simdata$X, disOptns = "correlation")
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_lt(ROE, 1e-2)
})

test_that("Stringing works for simulated example with moderate error using correlation metric",{
  set.seed(3)
  simdata = stringing_sim1(SNR = 4)
  stringingfit = Stringing(simdata$X, disOptns = "correlation")
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_lt(ROE, 1e-2)
})

test_that("Stringing works for simulated example with small error using euclidean distance",{
  set.seed(4)
  simdata = stringing_sim1(SNR = 10)
  stringingfit = Stringing(simdata$X, disOptns = "euclidean")
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_lt(ROE, 1e-2)
})

test_that("Stringing works for simulated example with small error using euclidean distance with standardization",{
  set.seed(4)
  simdata = stringing_sim1(SNR = 10)
  stringingfit = Stringing(simdata$X, disOptns = "euclidean", standardize = TRUE)
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_equal(ROE, 0)
})

test_that("Stringing works for simulated example with small error using spearman correlation",{
  set.seed(5)
  simdata = stringing_sim1(SNR = 10)
  stringingfit = Stringing(simdata$X, disOptns = "spearman")
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  # calculate Relative Order Error (ROE)
  ERp = (simdata$p-1)*(simdata$p+1)/3 # mean error for randomly sampled orders
  ROE = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))/ERp
  # the stringing function should recover the random ordering perfectly
  expect_lt(ROE, 1e-2)
})
