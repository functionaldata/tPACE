GenerateFunctionalData <-function(N, M, mu=NULL, lambda=NULL, k = 2,  basisType='cos'){
  
  if(N <2){
    stop("Sampes of size 1 are irrelevant.")
  }  
  if(M <20){
    stop("Dense samples with less than 20 observations per subject are irrelevant.")
  }
  s <- seq(0,1,length.out = M)
  if(is.null(mu)){
    mu = rep(0,M);
  }
  if(length(mu) != M){
    stop("Make sure that 'M' and the number of points over which 'mu' is evaluated is the same.")
  }
  if(is.null(lambda)){
    lambda = seq(k,1,-1)
  }
  if(k != length(lambda)){
    stop("Make sure you provide 'lambda's for all 'k' modes of variation.")
  }
  if( !(basisType %in% c('cos','sin','fourier'))){
    stop("Make sure you provide a valid parametric basis.")
  } 
  
  Ksi <- apply(matrix(rnorm(N*k), ncol=k), 2, scale) %*% diag(lambda, k)
  Phi <- CreateBasis(pts= s, type= basisType, K = k)
  
  yTrue <- t(matrix(rep(mu,N), nrow=M)) + Ksi %*% t(Phi) 
  return(list(Y = yTrue, Phi = Phi) )
}
 
