#' Create a sparse Functional Data sample for a Gaussian Process
#' 
#' Functional data sample of size n, sparsely sampled from a Gaussian process
#' 
#' @param n number of samples to generate.
#' @param rdist a sampler for generating the random design time points within [0, 1].
#' @param sparsity A vector of integers. The number of observation per sample is chosen to be one of the elements in sparsity with equal chance.
#' @param muFun a function that takes a vector input and output a vector of the corresponding mean (default: zero function).
#' @param K scalar specifying the number of basis to be used (default: 2).
#' @param lambda vector of size K specifying the variance of each components (default: rep(1,K)).
#' @param sigma The standard deviation of the Gaussian noise added to each observation points.
#' @param basisType string specifying the basis type used; possible options are: 'sin', 'cos' and 'fourier' (default: 'cos') (See code of 'CreateBasis' for implementation details.)
#' @param CovFun an alternative specification of the covariance structure.
#'
#' @return TODO
#' @export

MakeSparseGP <- function(n, rdist=runif, sparsity=2:9,
                         muFun=function(x) rep(0, length(x)), 
                         K = 2, lambda = rep(1, K), sigma=0,
                         basisType='cos', CovFun=NULL) {
   
  if(n < 2){
      stop("Samples of size 1 are irrelevant.")
  }  
  if(!is.function(rdist)){
      stop("'rdist' needs to be a function.")
  }
  if(!is.function(muFun)){
      stop("'muFun' needs to be a function.")
  }
  if (!is.numeric(sigma) || sigma < 0) {
    stop("'sigma' needs to be a nonnegative number")
  }
  if (!is.null(CovFun) && 
      (!missing(lambda) || !missing(basisType) || !missing(K))) {
    stop('Specify the covariance structure either with CovFun or with K, lambda, and basisType')
  }
  if (missing(K) && !missing(lambda) && is.null(CovFun)) {
    K <- length(lambda)
  }
  if(is.null(CovFun) && K != length(lambda)){
      stop("Make sure you provide 'lambda's for all 'K' modes of variation.")
  }
  # if( !(basisType %in% c('cos','sin','fourier'))) {
      # stop("Make sure you provide a valid parametric basis.")
  # } 
  if (length(sparsity) == 1) {
    sparsity <- c(sparsity, sparsity) # avoid scalar case for sample()
  }
   
  Ni <- sample(sparsity, n, replace=TRUE)
  if (is.null(CovFun)) {
    Ksi <- apply(matrix(rnorm(n*K), ncol=K), 2, scale) %*% 
      diag(sqrt(lambda), length(lambda))
     
    samp <- lapply(seq_len(n), function(i) {
      ni <- Ni[i]
      ti <- sort(rdist(ni))
      Phii <- CreateBasis(K, ti, basisType)
      yi <- muFun(ti) + as.numeric(tcrossprod(Ksi[i, ], Phii))

      list(ti = ti, yi=yi)
    })
  } else { # !is.null(CovFun)

    if (!requireNamespace('MASS', quietly=TRUE)) {
      stop('{MASS} is needed if CovFun is specified')
    }

    M <- 51
    pts <- seq(0, 1, length.out=M) # For generating the true curves

    samp <- lapply(seq_len(n), function(i) {
      ni <- Ni[i]
      ti <- sort(rdist(ni))
      ti <- c(ti, pts)
      yi <- MASS::mvrnorm(1, muFun(ti), CovFun(ti))
      
      list(ti = ti[seq_len(ni)], yi=yi[seq_len(ni)], yCurve=yi[-seq_len(ni)])
    })
  }

  Lt <- lapply(samp, `[[`, 'ti')
  Ly <- lapply(samp, `[[`, 'yi')
  
  
  if (sigma > 0) {
    LyTrue <- Ly
    Ly <- lapply(LyTrue, function(x) x + rnorm(length(x), sd=sigma))
  }

  if (is.null(CovFun)) {
    res <- list(Ly=Ly, Lt=Lt, xi=Ksi, Ni=Ni)
  } else {
    res <- list(Ly=Ly, Lt=Lt, yCurve=lapply(samp, `[[`, 'yCurve'), Ni=Ni)
  }

  if (sigma > 0) {
    res <- append(res, list(LyTrue=LyTrue))
  }
  
  return(res)
}
 
