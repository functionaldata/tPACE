devtools::load_all()
library(testthat)

# as in test_Rmullwlsk.R
try(silent=TRUE, load('data/InputFormMllwlskInCpp.RData'))
#try(silent=TRUE, load('tPACE/data/InputFormMllwlskInCpp.RData'))

IN = InputFormMllwlskInCpp
if(1==1){
  
  test_that('lwls2d interface is correct using xout1 and xout2', {
    AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), FALSE)
    BB = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='gauss',win=rep(1,38), FALSE)
    CC = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), FALSE)
    DD = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='epan',win=rep(1,38), FALSE)
    expect_equal(lwls2d(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid), AA)
    expect_equal(lwls2d(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn), BB)
    expect_equal(lwls2d(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid), CC)
    expect_equal(lwls2d(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn), DD)
  })
  
  
  test_that('lwls2d interface is correct using xout', {
    AA = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), FALSE))
    BB = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='gauss',win=rep(1,38), FALSE))
    CC = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), FALSE))
    DD = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='epan',win=rep(1,38), FALSE))
    expect_equal(lwls2d(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout=cbind(IN$regGrid, IN$regGrid)), AA, tolerance=0.05)
    expect_equal(lwls2d(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout=cbind(sort(unique(IN$tPairs[, 1])), sort(unique(IN$tPairs[, 2])))), BB)
    expect_equal(lwls2d(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout=cbind(IN$regGrid, IN$regGrid)), CC, tolerance=0.05)
    expect_equal(lwls2d(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout=cbind(sort(unique(IN$tPairs[, 1])), sort(unique(IN$tPairs[, 2])))), DD)
  })
  
}

if(1==1){
  
  test_that('lwls2d parallel works alright with a larger raw covariance',{
    
    set.seed(123)
    N = 200;
    M = 100;
    # Define the continuum
    s = seq(0,10,length.out = M)
    # Define the mean and 2 eigencomponents
    meanFunct <- function(s) s  + 10*exp(-(s-5)^2)
    eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
    eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
    
    # Create FPC scores
    Ksi = matrix(rnorm(N*2), ncol=2);
    Ksi = apply(Ksi, 2, scale)
    Ksi = Ksi %*% diag(c(5,2))
    
    # Create Y_true
    yTrue = Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))
    
    # Create sparse sample / Each subject has one to five readings (median: 3);
    ySparse = sparsify(yTrue, s, c(1:5))
    
    # Give your sample a bit of noise
    ySparse$yNoisy = lapply( ySparse$yList, function(x) x +  0.025*rnorm(length(x)))
    library(doParallel) 
    BB = GetRawCov(ySparse$yNoisy,ySparse$tList, s, meanFunct(s),'Sparse', TRUE)
    M = 9 
    BB$bw = c(3,4)
    ret  = lwls2d(c(3,4), kern="epan", xin =  (BB$tPairs), yin= (BB$cxxn), xout1= seq(1,9,length.out=M), xout2= seq(1,9,length.out=M) )
    
    parSmooth2 <- function(bw, xin, yin, kernel_type, win, nc, xout){ 
      if ( dim(xin)[1] != 2) {stop("You need to have xin input in a 2-by-n format")}
      if ( length(xout) < nc * 2 ){
        
      }
      N = length(xout) 
      breakPoints =   sort( N-  round(sqrt((1:(nc-1))/nc) * N))
      m2 <- Matrix::Matrix(0, nrow = N, N, sparse = TRUE)
      i_indx = c( breakPoints,  N ) # 1 3 7
      j_indx = c( 1, breakPoints+1) # 1 2 4 
      for (ij in 1:nc){ 
        q = i_indx[ij]
        p = j_indx[ij]
        m2[(p:q) , (p:N)] =  Rmullwlsk(bw, xin, cxxn= yin, xgrid=xout[p:N], ygrid=xout[p:q], 
                                       kernel_type=kernel_type, win=rep(1,length(yin)), FALSE, FALSE)
      } 
      m3 = as.matrix(m2) + t(as.matrix(m2))
      diag(m3) = 0.5 * diag(m3)
      return(m3)
    } 
     
    CCp  = parSmooth2(c(3,4), t(BB$tPairs), (BB$cxxn), 'epan', rep(1,length(BB$cxxn)), 3,  seq(1,9,length.out=M)) 
    CCp2 =   Rmullwlsk(c(3,4), kern="epan", tPairs=t(BB$tPairs), cxxn=(BB$cxxn), xgrid=seq(1,9,length.out=M), 
                       ygrid=seq(1,9,length.out=M), win=rep(1,1608),FALSE )  
    
    retP  = lwls2d(c(3,4), kern="epan", xin =  (BB$tPairs), yin= (BB$cxxn), xout1= seq(1,9,length.out=M), xout2= seq(1,9,length.out=M), userNumCores = 3 )

    expect_equivalent(CCp, ret)
    expect_equivalent(CCp, retP)
    expect_equivalent(CCp, CCp2) 
  })
  
}

if(1==1){
  
  test_that('lwls2d interface is correct using xout1 and xout for crosscovariances', {
    
    tPairs = matrix(c(1,3,4,1,2,1,1,2,3,4), ncol=2)
    
    AA = lwls2d(bw=c(0.5,0.5), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    AA2 = lwls2d(bw=c(0.5,0.5), kern ='gauss', xin=cbind(tPairs[, 2], tPairs[, 1]), yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    BB = lwls2d(bw=c(0.5,5.0), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    CC = lwls2d(bw=c(5.0,5.0), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    DD = lwls2d(bw=c(5.0,5.0), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,40)), crosscov=TRUE); 
    ZZ = lwls2d(bw=c(5.0,0.5), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4,4.5)), crosscov=TRUE); 
    # MATLAB equiv.
    # [invalid, AA]= mullwlsk_2([0.5,0.5], 'gauss', [1 3 4 1 2; 1 1 2 3 4], [1 2 3 4 5]', [1 1 1 1 1], [1  4], [1  4])
    # [invalid, BB]= mullwlsk_2([0.5,5.0], 'gauss', [1 3 4 1 2; 1 1 2 3 4], [1 2 3 4 5]', [1 1 1 1 1], [1  4], [1  4])
    # [invalid, CC]= mullwlsk_2([5.0,5.0], 'gauss', [1 3 4 1 2; 1 1 2 3 4], [1 2 3 4 5]', [1 1 1 1 1], [1  4], [1  4])
    # [invalid, DD]= mullwlsk_2([5.0,5.0], 'gauss', [1 3 4 1 2; 1 1 2 3 4], [1 2 3 4 5]', [1 1 1 1 1], [1  4], [1 40])
    # [invalid, ZZ]= mullwlsk_2([5.0,0.5], 'gauss', [1 3 4 1 2; 1 1 2 3 4], [1 2 3 4 5]', [1 1 1 1 1], [1  4], [1 4 4.5]); sum( ZZ(:))
    
    expect_equal(AA, t(AA2))
    expect_equal(13.997323601735092, sum(AA), tolerance=1e-9)
    expect_equal(13.498112821557918, sum(BB), tolerance=1e-9)
    expect_equal(13.669203283501956, sum(CC), tolerance=1e-9)
    expect_equal(89.458361008948557, sum(DD), tolerance=1e-9)
    expect_equal(24.498113242794656, sum(ZZ), tolerance=1e-9)
  })
  
}




