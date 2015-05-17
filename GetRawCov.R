GetRawCov <- function(y,t,out1new, mu, regular, error){
#  obtain raw covariance
#  Input y :       1*n cell array of the observed repeated measurements from n subjects
#  Input t :       1*n cell array of the observed time points from n subjects
#  Input out1new:  1*m vector of time points correspond to mu
#  Input mu:       1*m vector of fitted mean functions from Step I, corresponding to
#                 pooled unique time points from t
#  Input regular: Output of IsRegular()
#  Input error:    TRUE with measurement error assumption
#                  FALSE without measurement error assumption
# 
#  Output res: a list that contains tpairn, cxxn, indx,win and cyy
#     tpairn:  2 * N  vector denotes the  pairs of time points for subject 
#                 concatenating as two vectors 
#                if error = 1, all (t_ij, t_ij) will be removed  
#       cxxn:    1 * N vector of raw covariance corresponding to tpairn      
#       indx:    1 * N vector of indices for each subject
#        win:    1 * N weight matrix for the 2-D smoother for covariance function
#        cyy:    1 * M vector of raw covariance corresponding to all pairs of time points,
#                i.e., it is the same as cxxn if error = 0

 
  ncohort <- length(y);
  out1 <- sort(unique(unlist(t)))
  mu <- mapX1d(x = out1new, y = mu, newx = out1);
  count <- NULL
  indx = NULL 

  if(regular == 'Sparse'){
  
    Ys = lapply(X = y, FUN=meshgrid)
    Xs = lapply(X = t, FUN=meshgrid)

    xx1 = unlist(do.call(rbind, lapply(Xs, '[', 'X')) )
    xx2 = unlist(do.call(rbind, lapply(Xs, '[', 'Y')) ) 
    yy2 = unlist(do.call(rbind, lapply(Ys, '[', 'Y')) )
    yy1 = unlist(do.call(rbind, lapply(Ys, '[', 'X')) )

    id1 = apply(X= sapply(X=xx1, FUN='==',  ...=sort(unique(xx1)) ),MARGIN=2, FUN=which)
    id2 = apply(X= sapply(X=xx2, FUN='==',  ...=sort(unique(xx2)) ),MARGIN=2, FUN=which)
    cyy = ( yy1 - mu[ id1]) * (yy2 - mu[id2] )

    indx = unlist(sapply( 1:10, function(x) rep(x,  (unlist(lapply(length, X= y))[x])^2) ))

    tpairn = matrix( c(xx1, xx2), length(xx1),2);

    if(error){
      tneq = which(xx1 != xx2)
      indx = indx[tneq];
      tpairn = tpairn[tneq,];
      cxxn = cyy[tneq];     
    }else{
      cxxn = cyy;     
    }

    win = ones(1, length(cxxn));

  }else if(regular == 'Dense'){
    
    yy = t(matrix(unlist(y), length(y[[1]]), ncohort))
    MU = t(matrix( rep(mu, times=length(y)), ncol=length(y)))
    t1 = t[[1]]
    yy = yy - MU;
    cyy = t(yy) %*% yy / ncohort
    cyy = as.vector(t(cyy))
    cxxn = cyy;
    xxyy = meshgrid(t1);

    tpairn =  t(matrix( c(c(xxyy$X), c(xxyy$Y)), ncol = 2))

    if(error){
      tneq = which(tpairn[1,] != tpairn[2,])
      tpairn = tpairn[,tneq];
      cxxn = cyy[tneq];     
    }else{
      cxxn = cyy;     
    }

    win = ones(1, length(cxxn));

  }else if(regular == 'RegularWithMV'){
    stop("This is not implemented yet. Contact Pantelis!")
  }else {
    stop("Invalid 'regular' argument type")
  } 
    
  result <- list( 'tpairn'=tpairn, 'cxxn'=cxxn, 'indx'=indx, 'win'=win,'cyy'=cyy,'count'=count, 'error'=error, 'regular'=regular);  
 
  class(result) <- "RawCov"
  return(result)
}
