#' @title Confidence Intervals for Functional Linear Models.
#' @description Bootstrap pointwise confidence intervals for the coefficient functions in functional linear models.
#' @param Y Either an n-dimensional vector whose elements consist of scalar responses, or a list which contains functional responses in the form of a list LY and the time points LT at which they are observed (i.e., \code{list(Ly = LY,Lt = LT)}).
#' @param X A list of lists which contains the observed functional predictors list Lxj and the time points list Ltj at which they are observed. It needs to be of the form \code{list(list(Ly = Lx1,Lt = Lxt1),list(Ly = Lx2,Lt = Lxt2),...)}.
#' @param level A number taking values in [0,1] determining the confidence level. Default: 0.95.
#' @param R An integer holding the number of bootstrap replicates. Default: 999.
#' @param optnsListY A list of options control parameters for the response specified by \code{list(name=value)}. See 'Details' in FPCA.
#' @param optnsListX A list of options control parameters for the predictors specified by \code{list(name=value)}. See 'Details' in FPCA.

#' @details If measurement error is assumed, the diagonal elements of the raw covariance will be removed. This could result in highly unstable estimate 
#' if the design is very sparse, or strong seasonality presents. 
#' WARNING! For very sparse functional data, setting \code{measurementError=TRUE} is not recommended.
#' @return A list containing the following fields: 
#' \item{CI_alpha}{CI for the intercept function --- A data frame holding three variables: 
#' \code{CI_grid} --- the time grid where the CIs are evaluated,
#' \code{CI_lower} and \code{CI_upper} --- the lower and upper bounds of the CIs 
#' for the intercept function on \code{CIgrid}.}
#' \item{CI_beta}{ A list containing CIs for the slope functions --- the length of
#' the list is the same as the number of covariates. Each list contains the following fields:
#' A data frame holding three variables: \code{CI_grid} --- the time grid where the CIs are evaluated,
#' \code{CI_lower} and \code{CI_upper} --- the lower and upper bounds of the CIs 
#' for the coefficient function on \code{CIgrid} for \eqn{j = 1,2,\dots}.}
#' \item{level}{The confidence level of the CIs.}
#' @export


FLMCI <- function(Y, X, level = 0.95, R = 999, optnsListY = NULL, optnsListX = NULL){
  if (length(level) > 1) {
    level = level[1]
    warning("The input level has more than 1 element; only the first one is used.")
  }
  if (level < 0 | level > 1) {
    stop("Invalid input value of level.")
  }
  if (R %% 1 != 0 | R < 0) {
    stop("R should be a positive integer.")
  }
  
  if(is.list(Y)){
    n <- length(Y$Ly)
  }else{
    n <- length(Y)
  }
  
  p <- length(X)
  
  flm_est <- FLM(Y=Y, X=X, optnsListY=optnsListY, optnsListX=optnsListY)
  sparsityY <- flm_est$optnsListY$dataType
  sparsityX <- sapply(flm_est$optnsListX, function(x) x$dataType)
  if("Sparse" %in% c(sparsityY,sparsityX)){
    sparsity = "Sparse"
  }else{
    sparsity = "Dense"
  }
  beta_est <- flm_est$betaList
  workGridX <- flm_est$workGridX
  
  if(is.list(Y)){
    Y_temp <- range(unlist(Y$Lt))
    Y_l_ind <- which(sapply(Y$Lt, function(t){any(t %in% Y_temp[1])}))
    Y_u_ind <- which(sapply(Y$Lt, function(t){any(t %in% Y_temp[2])}))
  }
  
  X_temp <- lapply(
    X[sapply(X, is.list)], 
    function(v) {
      return(range(unlist(v[['Lt']])))
    }
  )
  X_l_ind <- lapply(
    1:length(X[sapply(X, is.list)]), 
    function(j) {
      return(which(sapply(X[[j]]$Lt, function(t){any(t %in% X_temp[[j]][1])})))
    }
  )
  X_u_ind <- lapply(
    1:length(X[sapply(X, is.list)]), 
    function(j) {
      return(which(sapply(X[[j]]$Lt, function(t){any(t %in% X_temp[[j]][2])})))
    }
  )
  
  # R times Boostrap
  betaMat <- lapply(1:R, function(b) {
    OutOfRange <- TRUE
    while (OutOfRange) {
      if(sparsity == "Sparse"){
        if(is.list(Y)){
          if(any(sapply(X, is.list))){
            ind <- c(sample(x = seq_len(n-2*(length(X[sapply(X, is.list)])+1)), n, replace = TRUE),
                     sample(x = Y_l_ind, 1),
                     sample(x = Y_u_ind, 1),
                     unlist(lapply(1:length(X[sapply(X, is.list)]), function(j){
                       c(sample(x = X_l_ind[[j]], 1),
                         sample(x = X_u_ind[[j]], 1))
                     }))
            )
          }else{
            ind <- c(sample(x = seq_len(n-2*(length(X[sapply(X, is.list)])+1)), n, replace = TRUE),
                     sample(x = Y_l_ind, 1),
                     sample(x = Y_u_ind, 1)
            )
          }
        }else{
          ind <- c(sample(x = seq_len(n-2*(length(X[sapply(X, is.list)])+1)), n, replace = TRUE),
                   unlist(lapply(1:length(X[sapply(X, is.list)]), function(j){
                     c(sample(x = X_l_ind[[j]], 1),
                       sample(x = X_u_ind[[j]], 1))
                   }))
          )
        }
        
        # check boundary of X
        temp_ind <- lapply(X[sapply(X, is.list)], function(v) {return(range(unlist(v[['Lt']][ind])))})
        cond <- sapply(1:length(X_temp), function(i){temp_ind[[i]][1] == X_temp[[i]][1] & temp_ind[[i]][2] == X_temp[[i]][2]})
        
        # check boundary of Y
        if(is.list(Y)){
          temp_ind <- range(unlist(Y$Lt[ind]))
          condY <- temp_ind[1] == Y_temp[1] & temp_ind[2] == Y_temp[2]
          cond <- c(cond, condY)
        }
        if(all(cond)){
          OutOfRange <- FALSE
        }
      }else{
        ind <- sample(x = seq_len(n), size = n, replace = TRUE)
        
        # check boundary of X
        temp_ind <- lapply(X[sapply(X, is.list)], function(v) {return(range(unlist(v[['Lt']][ind])))})
        cond <- sapply(1:length(X_temp), function(i){temp_ind[[i]][1] == X_temp[[i]][1] & temp_ind[[i]][2] == X_temp[[i]][2]})
        
        # check boundary of Y
        if(is.list(Y)){
          temp_ind <- range(unlist(Y$Lt[ind]))
          condY <- temp_ind[1] == Y_temp[1] & temp_ind[2] == Y_temp[2]
          cond <- c(cond, condY)
        }
        if(all(cond)){
          OutOfRange <- FALSE
        }
      }
    }
    if(is.list(Y)){
      Y_ind <- Y
      Y_ind$Lt <- Y$Lt[ind]
      Y_ind$Ly <- Y$Ly[ind]
    }else{
      Y_ind <- Y[ind]
    }
    
    X_ind <- X
    for(j in 1:p){
      if ( is.list(X_ind[[j]]) ) {
        X_ind[[j]]$Lt <- X[[j]]$Lt[ind]
        X_ind[[j]]$Ly <- X[[j]]$Ly[ind]
      }else{
        X_ind[[j]] <- X[[j]][ind]
      }
    }
    res <- FLM(Y = Y_ind, X = X_ind, optnsListY = optnsListY, optnsListX = optnsListX)
    length(res$betaList)
    if(is.list(Y)){
      return(list(alpha = res$alpha, beta = res$betaList, R2 = res$R2, workGridX = res$workGridX, workGridY = res$workGridY))
    }else{
      return(list(alpha = res$alpha, beta = res$betaList, R2 = res$R2, workGridX = res$workGridX))
    }
  })
  
  if(is.list(Y)){
    CI_alpha <- apply(t(sapply(1:R, function(b){
      betaMat[[b]]$alpha
    }, simplify = TRUE)), 2,
    stats::quantile, c((1-level)/2, 1-(1-level)/2))
    CI_alpha <- data.frame(Est = flm_est$alpha[1,], 
                           CI_lower = CI_alpha[1,], 
                           CI_upper = CI_alpha[2,], 
                           CI_grid = betaMat[[1]]$workGridY)
  }else{
    CI_alpha <- stats::quantile(sapply(1:R, function(b) betaMat[[b]]$alpha, simplify = TRUE), c((1-level)/2, 1-(1-level)/2))
    CI_alpha <- data.frame(Est = flm_est$alpha, 
                           CI_lower = CI_alpha[1], 
                           CI_upper = CI_alpha[2])
  }
  
  
  CI_beta <- lapply(1:p, function(j){
    if(is.list(Y)){
      if(is.list(X[[j]])){
        # beta matrix: row-workGridX and column-workGridY 
        ci_beta_df <-  data.frame( t(apply(t(sapply(1:R, function(b){
          as.vector(betaMat[[b]]$beta[[j]]) # matrix -> vector: arrange by column
        })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))))
        ci_beta <- list()
        ci_beta$Est <- flm_est$betaList[[j]]
        ci_beta$CI_lower <- matrix(ci_beta_df[,1], nrow = length(betaMat[[1]]$workGridX[[j]]), ncol = length(betaMat[[1]]$workGridY))
        ci_beta$CI_upper <- matrix(ci_beta_df[,2], nrow = length(betaMat[[1]]$workGridX[[j]]), ncol = length(betaMat[[1]]$workGridY))
        ci_beta$workGridX <- betaMat[[1]]$workGridX[[j]]
        ci_beta$workGridY <- betaMat[[1]]$workGridY
        return(ci_beta)
      }else{
        ci_beta_df <-  data.frame( t(apply(t(sapply(1:R, function(b){
          betaMat[[b]]$beta[[j]]
        })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))))
        names(ci_beta_df) <-  c("CI_lower", "CI_upper") 
        ci_beta_df$CI_grid <- betaMat[[1]]$workGridY
        ci_beta_df$Est <- flm_est$betaList[[j]][1,]
        return(ci_beta_df[,c(4,1,2,3)])
      }
    }else{
      if(is.list(X[[j]])){
        ci_beta_df <-  data.frame( t(apply(t(sapply(1:R, function(b){
          betaMat[[b]]$beta[[j]]
        })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))))
        names(ci_beta_df) <-  c("CI_lower", "CI_upper") 
        ci_beta_df$CI_grid <- betaMat[[1]]$workGridX[[j]]
        ci_beta_df$Est <- flm_est$betaList[[j]]
        return(ci_beta_df[,c(4,1,2,3)])
      }else{
        ci_beta_df <-  stats::quantile(sapply(1:R, function(b) betaMat[[b]]$beta[[j]]), c((1-level)/2, 1-(1-level)/2))
        names(ci_beta_df) <-  c("CI_lower", "CI_upper") 
        ci_beta_df <- c("Est"=flm_est$betaList[[j]], ci_beta_df)
        return(ci_beta_df)
      }
    }
  })
  names(CI_beta) <- sapply(1:p, function(j) { sprintf("CI_beta%d", j)})
  
  return(list(CI_alpha = CI_alpha, CI_beta = CI_beta, level = level))
  
}


