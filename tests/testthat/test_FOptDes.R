cat("\nTests for 'FOptDes'")
library(testthat)

test_that("optimal designs for trajectory recovery for dense data does not return any errors", { 
  set.seed(1)
  n <- 50
  pts <- seq(0, 1, by=0.05)
  sampWiener <- Wiener(n, pts)
  sampWiener <- MakeFPCAInputs(IDs = rep(c(1:n),length(pts)), tVec = rep(pts,each = n), yVec = as.vector(sampWiener))
  # global
  res <- FOptDes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, p=3,
                 isSequential=FALSE, RidgeCand = seq(0.1,1,0.1))
  # sequential optimization
  res <- FOptDes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, p=3,
                 isSequential=TRUE, RidgeCand = seq(0.1,1,0.1))
})

test_that("optimal designs for trajectory recovery for sparse data does not return any errors",{
  set.seed(1)
  n <- 50
  pts <- seq(0, 1, by=0.05)
  sampWiener <- Wiener(n, pts)
  sampWiener <- Sparsify(sampWiener, pts, 4:6)
  # global
  res <- FOptDes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, p=3,
                isSequential=FALSE, RidgeCand = seq(2,10,1))
  # sequential optimization
  resseq <- FOptDes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, p=3,
                 isSequential=TRUE, RidgeCand = seq(2,10,1))
})

test_that("optimal designs for response prediction for dense data does not return any errors", { 
  eifnMat <- function(K,t){
    l <- diff(range(t))
    a <- t[1]
    mat <- c()
    for(k in 1:K){
      mat <- cbind(mat,sqrt(2/l)*cos((k/l)*pi*(t-a)))
    }
    return(mat)
  }
  n=100
  reggrid = seq(0,1,0.05)
  mu=rep(0,length(reggrid))
  phi=eifnMat(K=10,t=reggrid)
  lambda=c(30, 20, 12, 8, 30/c(5:10)^2)
  K=10
  errorvar=0.25
  resp_errorvar = 0.25
  isDense = TRUE 
  Sparsity = 4:8
  DenseTrue <- matrix(rep(mu,n),byrow=TRUE,nrow=n,ncol=length(reggrid)) # n by length(RegGrid) matrix
  scores <- matrix(0,nrow=n,ncol=length(lambda)) # n by K matrix of FPC scores
  for(i in 1:ncol(scores)){ scores[,i] <- rnorm(n,mean=0,sd=sqrt(lambda[i])) }
  DenseTrue <- DenseTrue + t(phi %*% t(scores))
  # Generate independent measurement errors
  errorMat <- matrix(rnorm(n*length(reggrid),mean=0,sd=sqrt(errorvar)),nrow=n,ncol=length(reggrid))
  DenseObs <- DenseTrue + errorMat
  RespTrue <- scores %*% c(1,-2,1,-2,rep(0,K-4))
  RespObs <- c(RespTrue + rnorm(n, mean=0, sd=sqrt(resp_errorvar)))
  t <- list(); y <- list();
  for(i in 1:n){ t[[i]] <- reggrid; y[[i]] <- DenseObs[i,] }
  # global
  #res <- FOptDes(Ly=y, Lt=t, Resp=RespObs, p=3,
  #               isSequential=FALSE, RidgeCand = seq(0.1,1,0.1))
  # sequential optimization
  res <- FOptDes(Ly=y, Lt=t, Resp=RespObs, p=3,
                 isSequential=TRUE, RidgeCand = seq(0.1,1,0.1))
})

test_that("trajectory recovery case with p=1 for Brownian Motion on [0,1] returns theoretically correct optimal design", { 
  set.seed(1)
  n <- 100
  pts <- seq(0, 1, by=0.05)
  sampWiener <- Wiener(n, pts)
  fpcainput <- MakeFPCAInputs(IDs = rep(1:50,each=length(pts)), tVec = rep(pts, n), 
                              yVec = t(sampWiener))
  # global
  res <- FOptDes(Ly=fpcainput$Ly, Lt=fpcainput$Lt, p=1,
                 isSequential=FALSE, RidgeCand = 0.05)
  expect_equal(res$OptDes, 0.75)
})

test_that("medfly25 data example: optimal designs for response prediction for sparse data does not return error", { 
  data(medfly25)
  set.seed(1)
  medinput = MakeFPCAInputs(IDs = medfly25$ID, tVec = medfly25$Days, 
                            yVec = medfly25$nEggs)
  n = length(unique(medfly25$ID))
  respidx = seq(from = 1, by = 25, length.out = n)
  Resp = medfly25$remain_nEggs[respidx]
  res <- FOptDes(Ly = medinput$Ly, Lt = medinput$Lt, Resp = Resp, p = 2,
                 isSequential = FALSE, RidgeCand = seq(60,70,1))
})