#' Calculation of functional correlation between two simultaneously observed processes.
#'
#' @param x A list of function values corresponding to the first process.
#' @param y A list of function values corresponding to the second process.
#' @param Lt A list of time points for both \code{x} and \code{y}.
#' @param bw A numeric vector for bandwidth of length either 5 or 1, specifying the bandwidths for E(X), E(Y), var(X), var(Y), and cov(X, Y). If \code{bw} is a scalar then all five bandwidths are chosen to be the same. 
#' @param kern Smoothing kernel for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" (default: "gauss")
#' @param Tout Output time points. Default to the sorted unique time points. 
#'
#' @details \code{FCCor} calculate only the concurrent correlation corr(X(t), Y(t)) (note that the time points t are the same). It assumes no measurement error in the observed values.
#' @return A list with the following components:
#' \item{corr}{A vector of the correlation corr(X(t), Y(t)) evaluated at \code{Tout}.}
#' \item{Tout}{Same as the input Tout.}
#' \item{bw}{The bandwidths used for E(X), E(Y), var(X), var(Y), and cov(X, Y).}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' nGridIn <- 50
#' sparsity <- 1:5 # must have length > 1
#' bw <- 0.2
#' kern <- 'epan'
#' T <- matrix(seq(0.5, 1, length.out=nGridIn))
#' 
#' ## Corr(X(t), Y(t)) = 1/2
#' A <- Wiener(n, T)
#' B <- Wiener(n, T) 
#' C <- Wiener(n, T) + matrix((1:nGridIn) , n, nGridIn, byrow=TRUE)
#' X <- A + B
#' Y <- A + C
#' indEach <- lapply(1:n, function(x) sort(sample(nGridIn, sample(sparsity, 1))))
#' tAll <- lapply(1:n, function(i) T[indEach[[i]]])
#' Xsp <- lapply(1:n, function(i) X[i, indEach[[i]]])
#' Ysp <- lapply(1:n, function(i) Y[i, indEach[[i]]])
#' 
#' plot(T, FCCor(Xsp, Ysp, tAll, bw)[['corr']], ylim=c(-1, 1))
#' abline(h=0.5)
#' @export

FCCor <- function(x, y, Lt, bw=stop('bw missing'), kern='epan', Tout=sort(unique(unlist(Lt)))) {
  
  stopifnot(!is.null(x) && !is.null(y) && !is.null(Lt))
  
  Xvec <- unlist(x)
  Yvec <- unlist(y)
  tvec <- unlist(Lt)
  ord <- order(tvec)
  Xvec <- Xvec[ord]
  Yvec <- Yvec[ord]
  tvec <- tvec[ord]
  Tall <- sort(unique(unlist(Lt)))

  if (length(bw) != 1 && length(bw) != 5)
    stop('bw length incorrect.')
  if (is.numeric(bw) && length(bw) == 1)
    bw <- rep(bw, 5)
    
  muX <- Lwls1D(bw[1], kern, npoly=1L, nder=0L, xin=tvec, yin=Xvec, win=rep(1, length(tvec)), xout=Tall)
  muY <- Lwls1D(bw[2], kern, npoly=1L, nder=0L, xin=tvec, yin=Yvec, win=rep(1, length(tvec)), xout=Tall)
  names(muX) <- Tall
  names(muY) <- Tall
  
  Xcent <- Xvec - muX[as.character(tvec)]
  Ycent <- Yvec - muY[as.character(tvec)]
  varX <- Lwls1D(bw[3], kern, npoly=1L, nder=0L, xin=tvec, yin=Xcent^2, win=rep(1, length(tvec)), xout=Tout)
  varY <- Lwls1D(bw[4], kern, npoly=1L, nder=0L, xin=tvec, yin=Ycent^2, win=rep(1, length(tvec)), xout=Tout)
  covXY <- Lwls1D(bw[5], kern, npoly=1L, nder=0L, xin=tvec, yin=Xcent * Ycent, win=rep(1, length(tvec)), xout=Tout)

# The denominator variance esitmates may be negative sometimes. Set them to
  # NaN.
  varX[varX <= 0 ] <- NaN
  varY[varY <= 0 ] <- NaN

  if (any(is.nan(varX)) || any(is.nan(varY)))
    warning('NaN produced because the variance estimate is negative')

  res <- list(corr = covXY / sqrt(varX * varY), 
              Tout = Tout, 
              bw = bw)

  return(res)
}

