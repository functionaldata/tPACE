library(MASS)
library(testthat)
# devtools::load_all()

test_that('Regression works on margins', {
  set.seed(1000)
  
  n <- 100
  Lt <- list()
  Ly <- list()
  Lx1 <- list()
  Lx2 <- list()
  
  for (i in 1:n) {
    Ni <- sample(10:15,1)
    
    Lt[[i]] <- sort(runif(Ni,0,1))
    Lx1[[i]] <- runif(Ni,0,1)
    Lx2[[i]] <- runif(Ni,0,1)
    Ly[[i]] <- Lt[[i]]*(cos(2*pi*Lx1[[i]]) + sin(2*pi*Lx2[[i]])) + rnorm(Ni,0,0.1)
    
  }
  
  LLx <- list(Lx1,Lx2)
  
  gridT <- seq(0,1,length.out=41)
  x0 <- seq(0,1,length.out=51)
  x <- cbind(x0,x0)
  
  ht <- 0.1
  hx <- c(0.1,0.1)
  
  tvam <- TVAM(Lt,Ly,LLx,gridT=gridT,x=x,ht=ht,hx=hx,K='epan')
  
  tvam$tvamComp
  
  g0Sbf <- tvam$tvamMean
  gjSbf <- tvam$tvamComp
  
  # Check regression works on margins of time
  #time margin
  expect_equal(gjSbf[[1]][length(gridT),], cos(2*pi*x0), tolerance=0.2)
  expect_equal(gjSbf[[2]][length(gridT),], sin(2*pi*x0), tolerance=0.2)
  #x0 margin
  expect_equal(gjSbf[[1]][,1], gridT, tolerance=0.2)
  expect_equal(gjSbf[[1]][,length(x0)], gridT, tolerance=0.2)
})

test_that('Low noise is better medium noise is better than high noise', {
  changeNoise=function(sigma2){
    set.seed(1000)
    
    n <- 100
    Lt <- list()
    Ly <- list()
    Lx1 <- list()
    Lx2 <- list()
    
    for (i in 1:n) {
      Ni <- sample(10:15,1)
      
      Lt[[i]] <- sort(runif(Ni,0,1))
      Lx1[[i]] <- runif(Ni,0,1)
      Lx2[[i]] <- runif(Ni,0,1)
      Ly[[i]] <- Lt[[i]]*(cos(2*pi*Lx1[[i]]) + sin(2*pi*Lx2[[i]])) + rnorm(Ni,0,sigma2)
      
    }
    
    LLx <- list(Lx1,Lx2)
    
    gridT <- seq(0,1,length.out=41)
    x0 <- seq(0,1,length.out=51)
    x <- cbind(x0,x0)
    
    ht <- 0.1
    hx <- c(0.1,0.1)
    
    tvam <- TVAM(Lt,Ly,LLx,gridT=gridT,x=x,ht=ht,hx=hx,K='epan')
    return(tvam)
  }

  tvamVeryNoisy <- changeNoise(sigma2=0.5)
  tvamNoisy <- changeNoise(sigma2=0.2)
  tvamClean <- changeNoise(sigma2=0.01)
  
  gjSbfVeryNoisy <- tvamVeryNoisy$tvamComp
  gjSbfNoisy <- tvamNoisy$tvamComp
  gjSbfClean <- tvamClean$tvamComp
  
  #Comp 1
  expect_lt(sum((gjSbfClean[[1]][length(gridT),] - cos(2*pi*x0))^2),   sum((gjSbfNoisy[[1]][length(gridT),] - cos(2*pi*x0))^2)) 
  expect_lt(sum((gjSbfNoisy[[1]][length(gridT),] - cos(2*pi*x0))^2),   sum((gjSbfVeryNoisy[[1]][length(gridT),] - cos(2*pi*x0))^2)) 
  expect_lt(sum((gjSbfClean[[1]][length(gridT),] - cos(2*pi*x0))^2),   sum((gjSbfVeryNoisy[[1]][length(gridT),] - cos(2*pi*x0))^2)) 
  #Comp 2
  expect_lt(sum((gjSbfClean[[2]][length(gridT),] - sin(2*pi*x0))^2),   sum((gjSbfNoisy[[2]][length(gridT),] - sin(2*pi*x0))^2)) 
  expect_lt(sum((gjSbfNoisy[[2]][length(gridT),] - sin(2*pi*x0))^2),   sum((gjSbfVeryNoisy[[2]][length(gridT),] - sin(2*pi*x0))^2)) 
  expect_lt(sum((gjSbfClean[[2]][length(gridT),] - sin(2*pi*x0))^2),   sum((gjSbfVeryNoisy[[2]][length(gridT),] - sin(2*pi*x0))^2)) 
})

test_that('Regression surface is continuous', {
  set.seed(1000)
  
  n <- 100
  Lt <- list()
  Ly <- list()
  Lx1 <- list()
  Lx2 <- list()
  
  for (i in 1:n) {
    Ni <- sample(10:15,1)
    
    Lt[[i]] <- sort(runif(Ni,0,1))
    Lx1[[i]] <- runif(Ni,0,1)
    Lx2[[i]] <- runif(Ni,0,1)
    Ly[[i]] <- Lt[[i]]*(cos(2*pi*Lx1[[i]]) + sin(2*pi*Lx2[[i]])) + rnorm(Ni,0,0.1)
    
  }
  
  LLx <- list(Lx1,Lx2)
  
  gridT <- seq(0,1,length.out=41)
  x0 <- seq(0,1,length.out=51)
  x <- cbind(x0,x0)
  
  ht <- 0.1
  hx <- c(0.1,0.1)
  
  tvam <- TVAM(Lt,Ly,LLx,gridT=gridT,x=x,ht=ht,hx=hx,K='epan')
  
  tvam$tvamComp
  
  g0Sbf <- tvam$tvamMean
  gjSbf <- tvam$tvamComp
  
  #both components
  for(j in 1:2){
    #cts across time - jumps no larger than 0.2 across adjacent time points
    for(i in 1:(length(gridT)-1)){
     expect_equal(gjSbf[[j]][i,], gjSbf[[j]][i+1,], tol=0.2)
    }
    #cts across x0 - jumps no larger than 0.2 across adjacent x grid points
    for(i in 1:(length(x0)-1)){
      expect_equal(gjSbf[[j]][,i], gjSbf[[j]][,i+1], tol=0.2)
    }
  }
})

