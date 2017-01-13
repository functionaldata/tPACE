# Find optimal designs for trajectory recovery
BestDes_TR <- function(p, ridge, workGrid, Cov, isSequential=FALSE){
  # select optimal designs for trajectory recovery case, sequential method available
  if(isSequential == FALSE){ # Global Selection
    comblist <- utils::combn(1:length(workGrid),p)
    temps <- rep(0,ncol(comblist))
    for(i in 1:ncol(comblist)){  temps[i] <- TRCri(comblist[,i], ridge, Cov, workGrid)  }
    best <- sort(comblist[,min(which(temps==max(temps)))])
    return(list(best=best))
  } else { # Sequential optimization
    optdes <- c()
    for(iter in 1:p){
      candidx <- which(!((1:length(workGrid)) %in% optdes))
      seqcri <- rep(NA, length(candidx))
      for(i in 1:length(candidx)){
        tempdes <- sort(c(optdes,candidx[i]))
        seqcri[i] <- TRCri(tempdes, ridge, Cov, workGrid)
      }
      optdes <- sort(c(optdes, candidx[min(which(seqcri == max(seqcri)))]))
    }
    return(list(best=optdes,med=NULL)) # based on sequential selection
  }
}  

TRCri <- function(design, ridge, Cov, workGrid){
  # Optimization criterion for TR
  # Numerical Integration, equal to matrix multiplication if time grid is year
  design <- sort(design)
  RidgeCov <- Cov + diag(ridge, nrow(Cov))
  designcovinv <- solve(RidgeCov[design,design])
  if(length(design) > 1){
    trcri <- trapzRcpp(X=workGrid,Y=diag(t(Cov[design,])%*%designcovinv%*%(Cov[design,])))
  } else {
    trcri <- trapzRcpp(X=workGrid,Y=diag(Cov[design,]%*%designcovinv%*%(Cov[design,])))
  }
  return(trcri)
}
