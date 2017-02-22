# test CreateStringingPlot
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

test_that("CreateStringingPlot works",{
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
  CreateStringingPlot(stringingObj = stringingfit, subset = 1:10)
})

test_that("CreateStringingPlot works",{
  set.seed(1)
  simdata = stringing_sim1(SNR = Inf)
  stringingfit = Stringing(simdata$X, disOptns = "euclidean", standardize = TRUE)
  # check with simulated data to see if reversal of order is needed
  diff_norev = sum(abs(simdata$order[stringingfit$stringedPos] - 1:simdata$p))
  diff_rev = sum(abs(simdata$order[stringingfit$stringedPos] - simdata$p:1))
  if(diff_rev <= diff_norev){
    stringingfit$stringedPos = rev(stringingfit$stringedPos)
    stringingfit$Ly = lapply(stringingfit$Ly, rev)
  }
  CreateStringingPlot(stringingObj = stringingfit, subset = 1:10)
})
