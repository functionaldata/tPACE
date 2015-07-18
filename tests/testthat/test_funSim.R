 
  N = 333;
  set.seed(123)
  M = 100;
  
  # Define continuum
  s = seq(0,10,length.out = M)
 
  # Define mean and 2 eigencomponents
  meanFunct <- function(s) s  + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(c(5,2))
 
  # Create Y_true
  yTrue = Ksi %*% t( matrix(c( eigFunct1(s), eigFunct2(s)), ncol =2)) + t( matrix( rep(meanFunct(s),N), nrow= M))
  
  # Make a quick plot of what we have
  x11()
  par(mfrow=c(2,2))  
  plot(s,meanFunct(s), xlab='s',ylab='',  main="Mean Curve"); 
  plot(s,eigFunct1(s), xlab='s',ylab='', main="First Eigenfunc."); 
  plot(s,eigFunct2(s), xlab='s',ylab='',  main="Second Eigenfunc."); 
  matplot(t(yTrue), type='l', xlab='s',ylab='', main="True Sample Curves")
  
  # Create sparse sample
  ySparse16 = sparsify(yTrue, s, 1:16)
  ySparse08 = sparsify(yTrue, s, 1:8)
  ySparse04 = sparsify(yTrue, s, 1:4)
    
  # Give it a bit of noise
  ySparse16$yNoisy = lapply( ySparse16$yList, function(x) x +  0.5*rnorm(length(x)))
  ySparse04$yNoisy = lapply( ySparse04$yList, function(x) x +  0.5*rnorm(length(x)))
  ySparse08$yNoisy = lapply( ySparse08$yList, function(x) x +  0.5*rnorm(length(x)))
    
  A = FPCA(ySparse16$yNoisy, t= ySparse16$tList )
  B = FPCA(ySparse08$yNoisy, t= ySparse08$tList )
  C = FPCA(ySparse04$yNoisy, t= ySparse04$tList )
  QQ = makePACEinputs(IDs = rep(1:N, each=M),tVec=rep(s,N), t(yTrue) )
  D = FPCA(QQ$Ly, QQ$Lt) 

  x11()
  par(mfrow=c(2,4))   
  matplot((A$phi[,1:3]), type='l', xlab='s',ylab='', main='Eigenfunc. (median 9 p.)' )
  matplot((B$phi[,1:3]), type='l', xlab='s',ylab='', main='Eigenfunc. (median 5 p.)' )
  matplot((C$phi[,1:3]), type='l', xlab='s',ylab='', main='Eigenfunc. (median 3 p.)')
  matplot((D$phi[,1:2]), type='l', xlab='s',ylab='', main='Eigenfunc. (dense)')
  matplot((A$mu), type='l', xlab='s',ylab='', main='Mean (median 9 p.)' )
  matplot((B$mu), type='l', xlab='s',ylab='', main='Mean (median 5 p.)' )
  matplot((C$mu), type='l', xlab='s',ylab='', main='Mean (median 3 p.)' )
  matplot((D$mu), type='l', xlab='s',ylab='', main='Mean (dense)' )
  

  x11()
  par(mfrow=c(2,2))  
 createDesignPlot(t=ySparse16$tList, obsGrid= A$obsGrid, isColorPlot=FALSE, noDiagonal= TRUE,yname= '9 p.') 
 createDesignPlot(t=ySparse08$tList, obsGrid= B$obsGrid, isColorPlot=FALSE, noDiagonal= TRUE,yname= '5 p.') 
 createDesignPlot(t=ySparse04$tList, obsGrid= C$obsGrid, isColorPlot=FALSE, noDiagonal= TRUE,yname= '3 p.') 
 createDesignPlot(t=QQ$Lt, obsGrid= D$obsGrid, isColorPlot=FALSE, noDiagonal= FALSE,yname= 'dense') 
  
  
  
