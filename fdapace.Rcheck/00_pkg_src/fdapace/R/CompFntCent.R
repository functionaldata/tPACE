#####
##### centering component functions by using the marginal mean
#####

##### input variables: 
#####   f: evaluated values of component functions at estimation grid (N*d matrix)
#####   j: index of centering for the j-th component function (scalar)
#####   x: estimation grid (N*d matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (2-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   NW marginal regression function kernel estimators at estimation grid (N*d matrix)


# centering
CompFntCent <- function(f,j,x,MgnJntDens){
  
  fj <- f[,j]
  xj <- x[,j]
  
  pMatMgn <- MgnJntDens$pMatMgn
  
  tmp1 <- pMatMgn[,j]
  tmp <- fj-trapzRcpp(sort(xj),(fj*tmp1)[order(xj)])
  
  return(tmp)
  
}

