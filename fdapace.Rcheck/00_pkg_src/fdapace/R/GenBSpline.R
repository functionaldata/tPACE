# B-spline basis generator on [0,1] with equally spaced interior knots

GenBSpline <- function(x,nIntKnot=NULL,order=NULL) {
  
  # x: n-dimensional vector for evaluation
  # nIntKnot: a number of interior knots on (0,1) (scalar)
  # order: the order of B-spline basis function (scalar)
  
  if (is.null(nIntKnot)==TRUE) {
    nIntKnot <- 10
  }
  if (is.null(order)==TRUE) {
    order <- 3
  }
  if (nIntKnot < order) {
    stop('The number of knots should be greater than the order of B-spline basis.')
  }
  
  kOrder <- 0
  
  n <- length(x)
  t0 <- seq(0,1,length.out=(nIntKnot+2))[-c(1,nIntKnot+2)]
  
  newB <- matrix(0,nrow=n,ncol=(nIntKnot+1))
  t <- c(0,t0,1)
  
  for (i in 1:n) {
    for (k in 1:(nIntKnot)) {
      if (x[i] >= t[k] && x[i] < t[k+1]) {
        newB[i,k] <- 1
      }
      if (x[i]>=t[nIntKnot+1])
        newB[i,nIntKnot+1] <- 1
    }
  }
  
  if (order == 0) {
    
    return(B=newB)
    
  } else {
    
    while (kOrder < order) {
      
      kOrder <- kOrder + 1
      
      oldB <- cbind(rep(0,n),newB,rep(0,n))
      newB <- matrix(0,nrow=n,ncol=(nIntKnot+1+kOrder))
      t <- c(0,t,1)
      
      for (i in 1:n) {
        
        newB[i,1] <- (t[1+1+kOrder]-x[i])/(t[1+1+kOrder]-t[1+1])*oldB[i,1+1]
        
        for (k in 2:(nIntKnot+kOrder)) {
          newB[i,k] <- (x[i]-t[k])/(t[k+kOrder]-t[k])*oldB[i,k] + (t[k+1+kOrder]-x[i])/(t[k+1+kOrder]-t[k+1])*oldB[i,k+1]
        }
        
        newB[i,nIntKnot+1+kOrder] <- (x[i]-t[nIntKnot+1+kOrder])/(t[nIntKnot+1+2*kOrder]-t[nIntKnot+1+kOrder])*oldB[i,nIntKnot+1+kOrder]
      }
    }
    
    return(B=newB)
    
  }
}

# x <- sample(seq(0,1,length.out=201),201)
# 
# nIntKnot <- 10
# order <- 3
# 
# B <- GenBSpline(x,nIntKnot,order)
# # B <- GenBSpline(x)
# 
# plot(sort(x),B[order(x),1],type='l',col=1,ylim=c(min(B),max(B)))
# for (k in 2:ncol(B)) {
#   points(sort(x),B[order(x),k],type='l',col=k)
# }
# abline(v=seq(0,1,length.out=(nIntKnot+2)),col=8)
# B
# dim(B)