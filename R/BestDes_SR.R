# Find optimal designs for scalar response prediction
BestDes_SR <- function(p, ridge, workGrid, Cov, CCov, isSequential=FALSE){
  # select optimal designs for regression case, sequential method available
  if(isSequential == FALSE){
    comblist <- utils::combn(1:length(workGrid), p)
    temps <- rep(0,ncol(comblist))
    for(i in 1:ncol(comblist)){  temps[i] <- SRCri(comblist[,i], ridge, Cov, CCov)  }
    best <- sort(comblist[,min(which(temps==max(temps)))])
    return(list(best=best))
  } else{ # sequential selection
    optdes <- c()
    for(iter in 1:p){
      candidx <- which(!((1:length(workGrid)) %in% optdes))
      seqcri <- rep(NA, length(candidx))
      for(i in 1:length(candidx)){
        tempdes <- sort(c(optdes,candidx[i]))
        seqcri[i] <- SRCri(tempdes, ridge, Cov, CCov)
      }
      optdes <- sort(c(optdes, candidx[min(which(seqcri == max(seqcri)))]))
    }
    return(list(best=optdes,med=NULL))
  }
}

SRCri <- function(design,ridge,Cov,CCov){
  # Optimization criterion for SR
  design <- sort(design)
  ridgeCov <- Cov + diag(ridge,nrow(Cov))
  srcri <- t(CCov[design]) %*% solve(ridgeCov[design,design]) %*% CCov[design]
  return(srcri)
}
