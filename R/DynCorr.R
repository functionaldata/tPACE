#' @title Dynamical Correlation
#' @description Calculates the Dynamical Correlation for 2 paired dense regular functional data observed on the same grid.
#' @param x a n by m matrix where rows representing subjects and columns representing measurements, missings are allowed.
#' @param y a n by m matrix where rows representing subjects and columns representing measurements, missings are allowed.
#' @param t a length m vector of time points where x,y are observed.
#' @return A length n vector of individual dynamic correlations. The dynamic correlation can be obtained by taking average of this vector. 
#' @examples
#' set.seed(10)
#' n=200             # sample size
#' t=seq(0,1,length.out=100)       # length of data
#' mu_quad_x=8*t^2-4*t+5
#' mu_quad_y=8*t^2-12*t+6
#' fun=rbind(rep(1,length(t)),-t,t^2)
#' z1=matrix(0,n,3)
#' z1[,1]=rnorm(n,0,2)
#' z1[,2]=rnorm(n,0,16/3)
#' z1[,3]=rnorm(n,0,4)
#' x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))
#' for (i in 1:n){
#'   x1_quad_error[i,]=mu_quad_x+z1[i,]%*%fun+rnorm(length(t),0,0.01)
#'   y1_quad_error[i,]=mu_quad_y+2*z1[i,]%*%fun +rnorm(length(t),0,0.01)
#' }
#' dyn1_quad=DynCorr(x1_quad_error,y1_quad_error,t) 
#' @references
#' \cite{Dubin J A, Müller H G. Dynamical correlation for multivariate longitudinal data (2005).  
#' Journal of the American Statistical Association 100(471): 872-881.}
#' \cite{Liu S, Zhou Y, Palumbo R, Wang, J.L. (2016).  Dynamical correlation: A new method for quantifying synchrony with 
#' multivariate intensive longitudinal data. Psychological methods 21(3): 291.}
#' @export

DynCorr = function(x,y,t){
  if(dim(x)[1] != dim(y)[1] | dim(x)[2] != dim(y)[2]){
    stop("dimension of x and y does not match!")
  }
  if(dim(x)[2] != length(t)){
    stop("dimension of x,y does not match with t!")
  }
  na =  sum(is.na(x)+is.na(y))           

  if (na>0) {
    for (i in 1:dim(x)[1]){ 
      x[i,]=approx(t,x[i,],xout=t,rule=2)$y      
      y[i,]=approx(t,y[i,],xout=t,rule=2)$y
    }
  }

  temp1_x=temp1_y=matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  temp2_x=temp2_y=matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  M_x=M_y=z=numeric()

  for (i in 1:dim(x)[1]){
    aver_x=trapzRcpp(t,x[i,])/(tail(t,1)-head(t,1))
    aver_y=trapzRcpp(t,y[i,])/(tail(t,1)-head(t,1))
    temp1_x[i,]=x[i,]-aver_x
    temp1_y[i,]=y[i,]-aver_y
  }

  M_x=colMeans(temp1_x)
  M_y=colMeans(temp1_y) 
  for (i in 1:dim(x)[1]){
    temp2_x[i,]=(temp1_x[i,]-M_x)/sqrt(trapzRcpp(t,(temp1_x[i,]-M_x)^2)/(tail(t,1)-head(t,1)))
    temp2_y[i,]=(temp1_y[i,]-M_y)/sqrt(trapzRcpp(t,(temp1_y[i,]-M_y)^2)/(tail(t,1)-head(t,1)))
    z[i]=trapzRcpp(t,temp2_x[i,]*temp2_y[i,])/(tail(t,1)-head(t,1)) 
  }
  return(z)
}
