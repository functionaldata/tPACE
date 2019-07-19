GetMeanCurve=function(Ly, Lt, optns = list()){
  smcObj <- GetSmcObj(Ly, Lt, optns)
  mu <- smcObj$mu
  return(mu)
}

