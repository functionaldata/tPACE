#' @title Bootstrap test of Dynamic Correlation
#' @description Perform one sample (H0: Dynamic correlation = 0) or two sample (H0:Dynamic_correlation_1 = Dynamic_correlation_2) bootstrap test 
#' of H_0: Dynamical Correlation=0. 
#' @param x1 a n by m matrix where rows representing subjects and columns representing measurements, missings are allowed.
#' @param y1 a n by m matrix where rows representing subjects and columns representing measurements, missings are allowed.
#' @param t1 a vector of time points where x1,y1 are observed.
#' @param x2 (optional if missing will be one sample test) a n by m matrix where rows representing subjects and columns representing measurements, missings are allowed.
#' @param y2 (optional if missing will be one sample test) a n by m matrix where rows representing subjects and columns representing measurements, missings are allowed.
#' @param t2 (optional if missing will be one sample test) a vector of time points where x2,y2 are observed.
#' @param B number of bootstrap samples.
#' @return a list of the following 
#' \item{stats}{Test statistics.}
#' \item{pval}{p-value of the test.}
#' @examples
#' n=20             # sample size
#' t=seq(0,1,length.out=100)       # length of data
#' mu_quad_x=8*t^2-4*t+5
#' mu_quad_y=8*t^2-12*t+6
#' fun=rbind(rep(1,length(t)),-t,t^2)
#' z1=matrix(0,n,3)
#' z1[,1]=rnorm(n,0,2)
#' z1[,2]=rnorm(n,0,16/3)
#' z1[,3]=rnorm(n,0,4)   # covariance matrix of random effects
#' x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))
#' for (i in 1:n){
#'   x1_quad_error[i,]=mu_quad_x+z1[i,]%*%fun+rnorm(length(t),0,0.01)
#'   y1_quad_error[i,]=mu_quad_y+2*z1[i,]%*%fun +rnorm(length(t),0,0.01)
#' }
#' bt_DC=Dyn_test(x1_quad_error,y1_quad_error,t,B=500) # using B=500 for speed consideration
#'
#' @references
#' \cite{Dubin J A, Müller H G. (2005)  Dynamical correlation for multivariate longitudinal data.  
#' Journal of the American Statistical Association 100(471): 872-881.}
#'
#' \cite{Liu S, Zhou Y, Palumbo R, Wang, J.L. (2016).  Dynamical correlation: A new method for quantifying synchrony with multivariate intensive 
#' longitudinal data.  Psychological Methods 21(3): 291.}
#' @export

Dyn_test = function(x1,y1,t1,x2,y2,t2,B=1000){
  n=dim(x1)[1]
  if (missing(x2)) {                                           ### one-sample test ###
    na1 = sum(is.na(x1)+is.na(y1))
    if(sum(is.na(x1[,1])+is.na(x1[,ncol(x1)]))+ sum(is.na(y1[,1])+is.na(y1[,ncol(y1)])) > 0){
      warning('extrapolations may make results unreliable')
    }
    if (na1>0) {                                      ### impute missing values by linear interpolation ###
      for (i in 1:n){ 
        x1[i,]=approx(t1,x1[i,],xout=t1,rule=2)$y      
        y1[i,]=approx(t1,y1[i,],xout=t1,rule=2)$y
      }
    }
    dyncor_1=DynCorr(x1,y1,t1)                           ### observed DC ###
    obs_1_stud=mean(dyncor_1)*sqrt(n)/sd(dyncor_1)               ### observed standardized version ###
    boot_1_stud=numeric()
    boot_x1=boot_y1=matrix(0,nrow=n,ncol=length(t1))
    
    for(b in 1:B){                                               
      idx=sample(c(1:n),replace = TRUE)                         
      boot_x1=x1[idx,]                                          ### bootstrap replicates ###
      boot_y1=y1[idx,]
      boot_dyncor_1=DynCorr(boot_x1,boot_y1,t1)                  ### DC based on bootstrap samples ###
      boot_1_stud[b]=mean(boot_dyncor_1-dyncor_1)*sqrt(n)/sd(boot_dyncor_1)    ### bootstrap standardized version ###       
    }
    
    emp.stat=obs_1_stud                                          ### bootstrap test statistic ###
    emp.pval=length(boot_1_stud[boot_1_stud>abs(obs_1_stud) | boot_1_stud< -abs(obs_1_stud)])/B    ### p-value based on bootstrap null distribution ###
    outp=list(stats=emp.stat, pval=emp.pval) 
    return(outp)
  }  else {                                                    ### two-sample paired test (similar as above) ###
    na1 = sum(is.na(x1)+is.na(y1))
    na2 = sum(is.na(x2)+is.na(y2))                      #check missing
    if (na1 > 0) {                                      ### impute missing values by linear interpolation ###
      for (i in 1:n){ 
        x1[i,]=approx(t1,x1[i,],xout=t1,rule=2)$y      
        y1[i,]=approx(t1,y1[i,],xout=t1,rule=2)$y
      }
    }
    if (na2 > 0) {                                      ### impute missing values by linear interpolation ###
      for (i in 1:n){ 
        x2[i,]=approx(t1,x1[i,],xout=t1,rule=2)$y      
        y2[i,]=approx(t2,y2[i,],xout=t2,rule=2)$y
      }
    }
    dyncor_1=DynCorr(x1,y1,t1)
    dyncor_2=DynCorr(x2,y2,t2)
    obs_1_stud=mean(dyncor_1)*sqrt(n)/sd(dyncor_1)
    obs_2_stud=mean(dyncor_2)*sqrt(n)/sd(dyncor_2)
    obsdiff=dyncor_2-dyncor_1
    obsdiff_stud=mean(obsdiff)*sqrt(n)/sd(obsdiff)
    boot_1_stud=boot_2_stud=boot_diff_stud=numeric()
    
    boot_x1=boot_y1=matrix(0,nrow=n,ncol=length(t1))
    boot_x2=boot_y2=matrix(0,nrow=n,ncol=length(t2))
    
    for(b in 1:B){
      idx=sample(c(1:n),replace = TRUE)
      boot_x1=x1[idx,]
      boot_y1=y1[idx,]
      boot_x2=x2[idx,]
      boot_y2=y2[idx,]
      boot_dyncor_1=DynCorr(boot_x1,boot_y1,t1)
      boot_dyncor_2=DynCorr(boot_x2,boot_y2,t2)
      boot_1_stud[b]=mean(boot_dyncor_1-dyncor_1)*sqrt(n)/sd(boot_dyncor_1)
      boot_2_stud[b]=mean(boot_dyncor_2-dyncor_2)*sqrt(n)/sd(boot_dyncor_2)
      boot_diff_stud[b]=mean(boot_dyncor_2-boot_dyncor_1)*sqrt(n)/sd(boot_dyncor_2-boot_dyncor_1)
    }
    
    emp.stat=obsdiff_stud
    emp.pval=length(boot_diff_stud[boot_diff_stud>abs(obsdiff_stud) | boot_diff_stud< -abs(obsdiff_stud)])/B
    outp=list(stats=emp.stat, pval=emp.pval) 
    return(outp)
  } 
}
