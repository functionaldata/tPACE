#' Create a Dense Functional Data sample for a Gaussian process
#' 
#' For a Gaussian process, create a dense functional data sample of size n over a [0,1] support.
#' 
#' @param n number of samples to generate
#' @param M number of equidistant readings per sample (default: 100)
#' @param mu vector of size M specifying the mean (default: rep(0,M))
#' @param K scalar specifying the number of basis to be used (default: 2)
#' @param lambda vector of size K specifying the variance of each components (default: rep(1,K))
#' @param sigma The standard deviation of the Gaussian noise added to each observation points.
#' @param basisType string specifying the basis type used; possible options are: 'sin', 'cos' and 'fourier' (default: 'cos') (See code of 'CreateBasis' for implementation details.)
#'
#' @return Y: X(t_{j}), Yn: noisy observations
#' @export

MakeGPFunctionalData <-function(n, M = 100, mu=rep(0,M), K = 2, lambda = rep(1,K), sigma=0, basisType='cos'){
   
  if(n <2){
      stop("Samples of size 1 are irrelevant.")
  }  
  if(M <20){
      stop("Dense samples with less than 20 observations per subject are irrelevant.")
  }
  if (!is.numeric(sigma) || sigma < 0) {
    stop("'sigma' needs to be a nonnegative number")
  }
  s <- seq(0,1,length.out = M)

  if(length(mu) != M){
      stop("Make sure that 'M' and the number of points over which 'mu' is evaluated is the same.")
  }
  # if(is.null(lambda)){
      # lambda = seq(K,1,-1)
  # }
  if(K != length(lambda)){
      stop("Make sure you provide 'lambda's for all 'K' modes of variation.")
  }
  # if( !(basisType %in% c('cos','sin','fourier'))){
      # stop("Make sure you provide a valid parametric basis.")
  # } 
   
  Ksi <- apply(matrix(rnorm(n*K), ncol=K), 2, scale) %*% diag(sqrt(lambda))
  Phi <- CreateBasis(pts= s, type= basisType, K = K)
   
  yTrue <- t(matrix(rep(mu,n), nrow=M)) + Ksi %*% t(Phi) 
  
  res <- list(Y = yTrue, Phi = Phi, xi=Ksi, pts=s)
  
  if (sigma > 0) {
    yNoisy <- yTrue + rnorm(n * M, sd=sigma)
    res <- c(res, list(Yn = yNoisy))
  }
  
  return(res)
 }
 
