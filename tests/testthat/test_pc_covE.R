devtools::load_all()

p0 <- CreateOptions(dataType='Sparse', error=TRUE, kernel='epan')
p1 <- CreateOptions(dataType='Dense', error=TRUE, kernel='epan')

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)

sampSparse <- wiener(100, pts, sparsify=20)
p0 <- SetOptions(sampSparse$yList, sampSparse$tList, p0)
sampDense <- wiener(100, pts, sparsify=21)
p1 <- SetOptions(sampDense$yList, sampDense$tList, p1)
  
mu3 <- rep(0, length(pts))

## Sparse Case
rcovSparse <- GetRawCov(sampSparse$yList, sampSparse$tList, pts, mu3, p0$dataType, error=p0$error)
tmpSparse <- GetSmoothedCovarSurface(sampSparse$yList, sampSparse$tList, mu3, pts, regGrid,
                                     p0, useBins=FALSE)
pccovE_Sparse <- pc_covE(pts, seq(0,1,length.out=nrow(tmpSparse$smoothCov)), 
                         bw_userCov = tmpSparse$bwCov, rotationCut = c(0, 1), kernel = p0$kernel, rcov = rcovSparse)
pccovE_Sparse$sigma2

## Dense Case
rcovDense <- GetRawCov(sampDense$yList, sampDense$tList, pts, mu3, p1$dataType, error=p1$error)
tmpDense <- GetSmoothedCovarSurface(sampDense$yList, sampDense$tList, mu3, pts, regGrid,
                                     p1, useBins=FALSE)
pccovE_Dense <- pc_covE(pts, seq(0,1,length.out=nrow(tmpDense$smoothCov)), 
                         bw_userCov = tmpDense$bwCov, rotationCut = c(0, 1), kernel = p1$kernel, rcov = rcovDense)
pccovE_Dense$sigma2

## Tests added by Xiongtao
## Dense
set.seed(1)
pts <- seq(0, 1, by=0.01)
regGrid <- seq(0, 1, by=0.02)

n <- 10000
sampDense <- wiener(n, pts)
sampDense <- sampDense + rnorm(n * length(pts), sd=10)
sampDense <- sparsify(sampDense, pts, length(pts))
p1 <- SetOptions(sampDense$yList, sampDense$tList, list(dataType='Dense', error=TRUE, kernel='epan'))
mu3 <- rep(0, length(pts))

rcovDense <- GetRawCov(sampDense$yList, sampDense$tList, pts, mu3, p1$dataType, error=p1$error)
tmpDense <- GetSmoothedCovarSurface(sampDense$yList, sampDense$tList, mu3, pts, regGrid, p1, useBins=FALSE)


pccovE_Dense <- pc_covE(pts, seq(0,1,length.out=nrow(tmpDense$smoothCov)), bw_userCov = tmpDense$bwCov, rotationCut = c(0, 1), kernel = p1$kernel, rcov = rcovDense)
sqrt(pccovE_Dense$sigma2)

## Sparse
set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)
n <- 100
sampSparse <- wiener(n, pts)
sampSparse <- sampSparse + rnorm(n * length(pts), sd=0.4)
sampSparse <- sparsify(sampSparse, pts, 2:7)
p2 <- SetOptions(sampSparse$yList, sampSparse$tList, list(dataType='Sparse', error=TRUE, kernel='epan'))
mu3 <- rep(0, length(pts))

rcovSparse <- GetRawCov(sampSparse$yList, sampSparse$tList, pts, mu3, p2$dataType, error=p2$error)
tmpSparse <- GetSmoothedCovarSurface(sampSparse$yList, sampSparse$tList, mu3, pts, regGrid, p2, useBins=FALSE)

pccovE_Sparse <- pc_covE(pts, seq(0,1,length.out=nrow(tmpSparse$smoothCov)), bw_userCov = tmpSparse$bwCov, rotationCut = c(0, 1), kernel = p2$kernel, rcov = rcovSparse)
sqrt(pccovE_Sparse$sigma2)
