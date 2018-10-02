#'@title Dynamical Correlation
#'@description Calculate Dynamical Correlation for 2 paired dense regular functional data observed on the same grid.

#'@param x a n by m matrix where rows representing subjects and columns representing measurements
#'@param y a n by m matrix where rows representing subjects and columns representing measurements
#'@param t a length m vector of time points where x,y are observed.
#' @return A length m vector of individual dynamic correlations

DynCorr = function(x,y,t){
  x = t(x)
  y = t(y)
  na =  sum(is.na(x)+is.na(y))             

  if (na>0) {
    for (i in 1:dim(x)[2]){ 
      x[,i]=approx(t,x[,i],xout=t,rule=2)$y      
      y[,i]=approx(t,y[,i],xout=t,rule=2)$y
    }
  }

  temp1_x=temp1_y=matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  temp2_x=temp2_y=matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  M_x=M_y=z=numeric()

  for (i in 1:dim(x)[2]){
    aver_x=trapzRcpp(t,x[,i])/(tail(t,1)-head(t,1))
    aver_y=trapzRcpp(t,y[,i])/(tail(t,1)-head(t,1))
    temp1_x[,i]=x[,i]-aver_x
    temp1_y[,i]=y[,i]-aver_y
  }

  M_x=rowMeans(temp1_x)
  M_y=rowMeans(temp1_y) 
  for (i in 1:dim(x)[2]){
    temp2_x[,i]=(temp1_x[,i]-M_x)/sqrt(trapzRcpp(t,(temp1_x[,i]-M_x)^2)/(tail(t,1)-head(t,1)))
    temp2_y[,i]=(temp1_y[,i]-M_y)/sqrt(trapzRcpp(t,(temp1_y[,i]-M_y)^2)/(tail(t,1)-head(t,1)))
    z[i]=trapzRcpp(t,temp2_x[,i]*temp2_y[,i])/(tail(t,1)-head(t,1)) 
  }
  return(z)
}