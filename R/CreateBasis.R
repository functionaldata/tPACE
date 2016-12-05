# Create an orthogonal basis of K functions in [0, 1], with nGrid points.
# type: a string in 'sin', 'cos', 'fourier', 'legendre01'
# Output: a K by nGrid matrix, each column containing an basis function.
CreateBasis <- function(K, pts=seq(0, 1, length.out=50), type='sin') {

  nGrid <- length(pts)
  possibleTypes <- c('cos', 'sin', 'fourier', 'legendre01', 'poly', 'unknown')
  type <- possibleTypes[pmatch(type, possibleTypes, nomatch=length(possibleTypes))]
  
  stopifnot(is.numeric(K) && length(K) == 1 && K > 0)
  
  if (type == 'cos') {
    res <- sapply(seq_len(K), function(k) 
      if (k == 1) {
        rep(1, nGrid)
      } else {
        sqrt(2) * cos((k - 1) * pi * pts)
      }
    )
  } else if (type == 'sin') {
    res <- sapply(seq_len(K), function(k) sqrt(2) * sin(k * pi * pts))
  } else if (type == 'fourier') {
    res <- sapply(seq_len(K), function(k) 
      if (k == 1) {
        rep(1, nGrid)
      } else if (k %% 2 == 0) {
        sqrt(2) * sin(k * pi * pts)
      } else {
        sqrt(2) * cos((k - 1) * pi * pts)
      }
    )
  } else if (type == 'legendre01') {
    # coefMat <- matrix(0, K, K)
    if (K == 1) {
      res <- matrix(1, length(pts), 1)
    } else if (K > 1) {
      coefMat <- sapply(seq_len(K), function(n) {
        coef <- rep(0, K)
        # # coef[1] <- (-1)^(n - 1)
        for (k in seq_len(n)) {
          coef[k] <- (-1)^(n - k) * choose(n - 1, k - 1) * 
                                    choose(n + k - 2, k - 1)
        }
        coef * sqrt(2 * n - 1)
      })
      xMat <- cbind(1, stats::poly(pts, K - 1, raw=TRUE))
      res <- xMat %*% coefMat
      # browser()
    }
  } else if (type == 'poly') {
    if (K == 1) {
      res <- matrix(1, length(pts), 1)
    } else if (K > 1) {
     res <- cbind(1, stats::poly(pts, K - 1, raw=TRUE))
    }
  } else if (type == 'unknown') {
    stop('unknown basis type')
  }

  res <- matrix(res, ncol=K) # prevent single length pts
  res
}
