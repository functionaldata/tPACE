#' Makke Gaussian Process Dense Functional Data sample                                                             
#' 
#' Make a Gaussian process dense functional data sample of size N over a [0,1] support.
#' 
#' @param N number of samples to generate
#' @param M number of equidistant readings per sample (default: 100)
#' @param mu vector of size M specifying the mean (default: rep(0,M))
#' @param k scalar specifying the number of basis to be used (default: 2)
#' @param lambda vector of size K specifying the variance of each components (default: rep(1,k))
#' @param basisType string specifiying the basis type used; possible options are: 'sin', 'cos' and 'fourier' (default: 'cos') (See code of 'CreateBasis' for implementation details.)
#'

MakeGPFunctionalData <-function(N, M = 100, mu=rep(0,M), k = 2, lambda = rep(1,k),  basisType='cos'){
   
  if(N <2){
      stop("Samples of size 1 are irrelevant.")
  }  
  if(M <20){
      stop("Dense samples with less than 20 observations per subject are irrelevant.")
  }
  s <- seq(0,1,length.out = M)

  if(length(mu) != M){
      stop("Make sure that 'M' and the number of points over which 'mu' is evaluated is the same.")
  }
  # if(is.null(lambda)){
      # lambda = seq(k,1,-1)
  # }
  if(k != length(lambda)){
      stop("Make sure you provide 'lambda's for all 'k' modes of variation.")
  }
  if( !(basisType %in% c('cos','sin','fourier'))){
      stop("Make sure you provide a valid parametric basis.")
  } 
   
  Ksi <- apply(matrix(rnorm(N*k), ncol=k), 2, scale) %*% diag(sqrt(lambda))
  Phi <- CreateBasis(pts= s, type= basisType, K = k)
   
  yTrue <- t(matrix(rep(mu,N), nrow=M)) + Ksi %*% t(Phi) 
  return(list(Y = yTrue, Phi = Phi, xi=Ksi, pts=s) )
 }
 
