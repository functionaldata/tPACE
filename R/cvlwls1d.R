cvlwls1d <- function(yy, tt, kernel, npoly, nder, dataType ){

  ncohort = length(t);
  tt  = unlist(t);
  xx  = unlist(yy);
  ind = unlist(lapply( 1:length(t), function(j) rep(j, times=length(t[[j]]))));

  ttn=tt;
  xxn=xx;
  a0=min(tt);
  b0=max(tt);
  rang = b0-a0;
  dstar = minb(tt, npoly+2);

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
        out=tt[ind==i];
        obs=xx[ind==i];
        win=rep(1,length(tt));
        win[ind==i]=0;        
        
        if (dataType=='Dense') {
            xxn=(ave*ncohort-t[[i]])/(ncohort-1);
            ttn=t[[1]];
            win=ones(1,length(t[[1]]));
        }

        mu= lwls1d(bw= bw[j], kern=kernel, npoly=npoly, nder= nder, xin = ttn, yin= xxn, xout=out)       

        # if invalid==0 {
        cv[j]=cv[j]+t(obs-mu)%*%(obs-mu);
        count[j]=count[j]+1;
        #}
    }
  }
  cv = cv[(count/ncohort>0.90)];
  bw = bw[(count/ncohort>0.90)];
  bopt = bw[(cv==min(cv))];
  
  return(bopt)

}

