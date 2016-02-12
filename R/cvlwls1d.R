CVLwls1D <- function(yy, t, kernel, npoly, nder, dataType ){

  # If 'yy' and 't' are vectors "cheat" and break them in 
  # a list of 10 elements
  if ( is.vector(yy) && is.vector(t) && !is.list(t) && !is.list(yy) ){
    if (length(t) < 21) {
      stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    }
    myPartition =   c(1:10, sample(10, length(t)-10, replace=TRUE));
    yy = split(yy, myPartition)
    t = split(t, myPartition)
    dataType = 'Sparse';
  }

  ncohort = length(t);
  tt  = unlist(t);
  xx  = unlist(yy);
  ind = unlist(lapply( 1:ncohort, function(j) rep(j, times=length(t[[j]]))));

  xxn = xx[order(tt)];
  ind = ind[order(tt)];
  ttn = sort(tt);
  
  a0=ttn[1];
  b0=ttn[length(ttn)];
  rang = b0-a0;
  dstar = Minb(tt, npoly+2);

  if (dataType != 'Dense'){
   h0 = 2.5*dstar;
  } else {
   h0 = dstar;
  }
  
  if (h0 > rang/4){
    h0 = h0*.75;
    warning(sprintf("Warning: the min bandwith choice is too big, reduce to %f !", (h0)  ))
  }    
  
  nbw = 11;
  bw = rep(0,nbw-1);
  n = length(unique(tt));
  
  for (i in 1:(nbw-1)){
    bw[i]=2.5*rang/n*(n/2.5)^((i-1)/(nbw-1)); # Straight from MATLAB
  }
  bw = bw-min(bw)+h0;
  
  ave = rep(0, length(t[[1]]));
  
  if (dataType == 'Dense'){
    for (i in 1:ncohort){
      ave = ave + t[[i]]/ncohort;
    }
  }

  cv = c();
  count = c();

  for (j in 1:(nbw-1)){
    cv[j]=0;
    count[j]=0;
    for (i in 1:ncohort){
      out=ttn[ind==i];
      obs=xxn[ind==i];
      win=rep(1,length(tt));
      win[ind==i] = NA;        
        
      if (dataType=='Dense') {
        xxn=(ave*ncohort-t[[i]])/(ncohort-1);
        ttn=t[[1]];
        win=pracma::ones(1,length(t[[1]]));    
        xxn = xxn[order(ttn)]
        ttn = sort(ttn)           
      }  
 
      mu = Lwls1D(bw= bw[j], kernel_type = kernel, npoly=npoly, nder= nder, xin = ttn, yin= xxn, xout=out, win = win)
      cv[j]=cv[j]+t(obs-mu)%*%(obs-mu);
      count[j]=count[j]+1;
    }
  }
  cv = cv[(count/ncohort>0.90)];
  bw = bw[(count/ncohort>0.90)];
  bopt = bw[(cv==min(cv))];
  
  return(bopt)

}

