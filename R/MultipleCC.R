
MultipleCC = function( expVarFunc, bwFunct = NULL, kernelType = 'gauss' ){
  
  Q =  length(expVarFunc)
  L =  length(expVarFunc[[1]]$workGrid)

  # Construct the matrix CCXX auto- and cross- covariances for the functional predictors
  CCXX = matrix( rep(0, ( L* Q)^2 ), nrow = L* Q ) 
  for(i in 2:( Q -0)){
    if( !all.equal(expVarFunc[[i]]$workGrid, expVarFunc[[1]]$workGrid) ){
      stop('You need the same workGrid for the functional variables') 
    }
  }
  for( i in 1:Q){
    for( j in i:Q){
      if(i == j){
    CCXX[ ((i-1)*L) + (1:L), ((j-1)*L) + (1:L) ] = expVarFunc[[j]]$smoothedCov
      } else {
    tempCCXX =  CrCovYX(Ly1 = expVarFunc[[i]]$inputData$y, Ly2 = expVarFunc[[j]]$inputData$y,
           Lt1 = expVarFunc[[i]]$inputData$t, Lt2 = expVarFunc[[j]]$inputData$t, kernelType = kernelType,
           Ymu1= expVarFunc[[i]]$mu,  Ymu2= expVarFunc[[j]]$mu, bw1 = bwFunct[1], bw2 = bwFunct[2]);
    CCXX[ ((i-1)*L) + (1:L), ((j-1)*L) + (1:L) ] = tempCCXX$smoothedCC;
    # And the transpose
    CCXX[ ((j-1)*L) + (1:L), ((i-1)*L) + (1:L) ] = t(tempCCXX$smoothedCC);
      }
    }  
  }
  return(CCXX)                
}
