#' Estimation of measurement error \code{sigma2} via Z Lin and JL Wang (2021+)
#'
#' @param y A list of \emph{n} vectors containing the observed values for each individual.
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' 
#' @return 
#' \item{sigma2}{Variance for measurement error.}
#' 
#' @references
#' \cite{Z Lin and JL Wang (2021+). "Mean and covariance estimation for functional snippets." JASA.}
#' @export

GetLWSigma2 <- function(y,t){
  if(is.list(y)) # irregular data
    {
      n <- length(t)
      select.sig2.bw <- function(Lt,Ly){#bandwidth selection
        n <- length(Lt)
        t.min <- min(unlist(Lt))
        t.max <- max(unlist(Lt))
        delta <- max(sapply(Lt,function(ts){max(ts)-min(ts)}))
        m <- mean(sapply(Lt,function(ts){length(ts)}))
        M <- n * m^2
        
        #res <- FPCA(Ly,Lt,optns=list(methodRho='vanilla',dataType=dataType,methodXi='CE'))#FPCA(Ly, Lt, list(usergrid=TRUE,error=FALSE,verbose=TRUE))#FPCA(Ly, Lt)
        #grid <- res$workGrid
        #mu <- res$mu
        rs<- GetMeanCurve(Ly,Lt)
        grid<- rs$workGrid
        mu <- rs$mu
        
        mu_est <- function(x,grid,mu){
          res <- c()
          for (i in 1:length(x)) {
            res[i] <- mu[which(min(abs(x[i]-grid))==abs(x[i]-grid))]
          }
          return(res)
        }
        mu.hat <- lapply(Lt, mu_est,grid=grid,mu=mu)
        
        tmp <- lapply(1:length(Lt),function(i){
          rr <- Ly[[i]] - mu.hat[[i]]
          rr^2
        })
        vn <- sqrt(mean(unlist(tmp)))
        
        h <- 0.29 * delta * vn * (M^(-1/5))
        
        max.it <- 1000
        it <- 0
        while(it < max.it) # h0 two small
          {
            it <- it + 1
            
            cnt <- sapply(1:n,function(i){
              tobs <- Lt[[i]]
              y <- Ly[[i]]
              m <- length(tobs)
              v1 <- 0
              
              if(m < 2) return(0)
              for(j in 1:m)
              {
                for(k in 1:m)
                {
                  if( (k!=j) && abs(tobs[j]-tobs[k])<h )
                  {
                    v1 <- v1 + 1
                  }
                }
              }
              return(v1)
            })
            
            cnt <- sum(cnt)
            if(cnt >= min(50,0.1*n*m*(m-1))) break
            else
            {
              h <- h*1.01
            }
        }
        return(h)
      }
      h <- select.sig2.bw(t,y)
      AB <- sapply(1:n,function(i){
        tobs <- t[[i]]
        y <- y[[i]]
        m <- length(tobs)
        v1 <- 0
        v2 <- 0
        
        if(m < 2) return(c(v1,v2))
        for(j in 1:m)
        {
          for(k in 1:m)
          {
            if( (k!=j) && abs(tobs[j]-tobs[k])< h )
            {
              v1 <- v1 + (y[j]-y[k])^2 / 2
              v2 <- v2 + 1
            }
          }
        }
        if(v2 == 0) return(0)
        else return(v1/v2)
      })
      AB <- unlist(AB)
      nz <- sum(AB!=0)
      if(nz==0) sig2 <- 0
      else sig2 <- sum(AB)/nz
  }else if(is.matrix(y)) # regular design
    {
      p <- ncol(y)
      n <- nrow(y)
      
      L1 <- y[,1:(p-2)]
      L2 <- y[,2:(p-1)]
      L3 <- y[,3:p]
      
      L <- L1+L3-2*L2
      L <- L^2
      sig2 <- mean(L)/(2*3)
  }else stop('unsupported data type')
  return(sig2)
}
