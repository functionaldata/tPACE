# devtools::load_all()
library(testthat)

# as in test_Rmullwlsk.R
try(silent=TRUE, load(system.file('testdata', 'InputFormMllwlskInCpp.RData', package='fdapace')))
#try(silent=TRUE, load(system.file('testdata', 'InputFormMllwlskInCpp.RData', package='fdapace')))

IN = InputFormMllwlskInCpp
if(1==1){
  
  test_that('Lwls2D interface is correct using xout1 and xout2', {
    AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), FALSE)
    BB = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='gauss',win=rep(1,38), FALSE)
    CC = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), FALSE)
    DD = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='epan',win=rep(1,38), FALSE)
    expect_equal(Lwls2D(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid), AA)
    expect_equal(Lwls2D(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn), BB)
    expect_equal(Lwls2D(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid), CC)
    expect_equal(Lwls2D(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn), DD)
  })
  
  
  test_that('Lwls2D interface is correct using xout', {
    AA = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), FALSE))
    BB = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='gauss',win=rep(1,38), FALSE))
    CC = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), FALSE))
    DD = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='epan',win=rep(1,38), FALSE))
    expect_equal(Lwls2D(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout=cbind(IN$regGrid, IN$regGrid)), AA, tolerance=0.05)
    expect_equal(Lwls2D(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout=cbind(sort(unique(IN$tPairs[, 1])), sort(unique(IN$tPairs[, 2])))), BB)
    expect_equal(Lwls2D(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout=cbind(IN$regGrid, IN$regGrid)), CC, tolerance=0.05)
    expect_equal(Lwls2D(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout=cbind(sort(unique(IN$tPairs[, 1])), sort(unique(IN$tPairs[, 2])))), DD)
  })
  
}

if(1==1){
  
  test_that('Lwls2D interface is correct using xout1 and xout for crosscovariances', {
    
    tPairs = matrix(c(1,3,4,1,2,1,1,2,3,4), ncol=2)
    
    AA = Lwls2D(bw=c(0.5,0.5), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    AA2 = Lwls2D(bw=c(0.5,0.5), kern ='gauss', xin=cbind(tPairs[, 2], tPairs[, 1]), yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    BB = Lwls2D(bw=c(0.5,5.0), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    CC = Lwls2D(bw=c(5.0,5.0), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4)), crosscov=TRUE);
    DD = Lwls2D(bw=c(5.0,5.0), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,40)), crosscov=TRUE); 
    ZZ = Lwls2D(bw=c(5.0,0.5), kern ='gauss', xin=tPairs, yin=c(1,2,3,4,5), xout1=as.numeric(c(1,4)), xout2=as.numeric(c(1,4,4.5)), crosscov=TRUE); 
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




