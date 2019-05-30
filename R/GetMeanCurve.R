GetMeanCurve=function(Ly, Lt, optns = list()){
  smcObj <- GetSmcObj(Ly, Ly, optns)
  mu <- smcObj$mu
  return(mu)
}
