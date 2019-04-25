# devtools::load_all()

p0 <- list(dataType='Sparse', error=TRUE, kernel='epan')
p1 <- list(dataType='Dense', error=TRUE, kernel='epan')

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)

sampSparse <- Wiener(100, pts, sparsify=20)
p0 <- SetOptions(sampSparse$Ly, sampSparse$Lt, p0)
sampDense <- Wiener(100, pts, sparsify=21)
p1 <- SetOptions(sampDense$Ly, sampDense$Lt, p1)
  
mu3 <- rep(0, length(pts))

## Sparse Case
rcovSparse <- GetRawCov(sampSparse$Ly, sampSparse$Lt, pts, mu3, p0$dataType, error=p0$error)
tmpSparse <- GetSmoothedCovarSurface(sampSparse$Ly, sampSparse$Lt, mu3, pts, regGrid,
                                     p0, useBinnedCov=FALSE)
pccovE_Sparse <- PC_CovE(pts, seq(0,1,length.out=nrow(tmpSparse$smoothCov)), 
                         bw_userCov = tmpSparse$bwCov, rotationCut = c(0, 1), kernel = p0$kernel, rcov = rcovSparse)
pccovE_Sparse$sigma2

## Dense Case
rcovDense <- GetRawCov(sampDense$Ly, sampDense$Lt, pts, mu3, p1$dataType, error=p1$error)
tmpDense <- GetSmoothedCovarSurface(sampDense$Ly, sampDense$Lt, mu3, pts, regGrid,
                                     p1, useBinnedCov=FALSE)
pccovE_Dense <- PC_CovE(pts, seq(0,1,length.out=nrow(tmpDense$smoothCov)), 
                         bw_userCov = tmpDense$bwCov, rotationCut = c(0, 1), kernel = p1$kernel, rcov = rcovDense)
pccovE_Dense$sigma2

## Tests added by Xiongtao
## Dense
set.seed(1)
pts <- seq(0, 1, by=0.01)
regGrid <- seq(0, 1, by=0.02)

n <- 10000
sampDense <- Wiener(n, pts)
sampDense <- sampDense + rnorm(n * length(pts), sd=10)
sampDense <- Sparsify(sampDense, pts, length(pts))
p1 <- SetOptions(sampDense$Ly, sampDense$Lt, list(dataType='Dense', error=TRUE, kernel='epan'))
mu3 <- rep(0, length(pts))

rcovDense <- GetRawCov(sampDense$Ly, sampDense$Lt, pts, mu3, p1$dataType, error=p1$error)
tmpDense <- GetSmoothedCovarSurface(sampDense$Ly, sampDense$Lt, mu3, pts, regGrid, p1, useBinnedCov=FALSE)


pccovE_Dense <- PC_CovE(pts, seq(0,1,length.out=nrow(tmpDense$smoothCov)), bw_userCov = tmpDense$bwCov, rotationCut = c(0, 1), kernel = p1$kernel, rcov = rcovDense)
sqrt(pccovE_Dense$sigma2)

## Sparse
set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)
n <- 100
sampSparse <- Wiener(n, pts)
sampSparse <- sampSparse + rnorm(n * length(pts), sd=0.4)
sampSparse <- Sparsify(sampSparse, pts, 2:7)
p2 <- SetOptions(sampSparse$Ly, sampSparse$Lt, list(dataType='Sparse', error=TRUE, kernel='epan'))
mu3 <- rep(0, length(pts))

rcovSparse <- GetRawCov(sampSparse$Ly, sampSparse$Lt, pts, mu3, p2$dataType, error=p2$error)
tmpSparse <- GetSmoothedCovarSurface(sampSparse$Ly, sampSparse$Lt, mu3, pts, regGrid, p2, useBinnedCov=FALSE)

pccovE_Sparse <- PC_CovE(pts, seq(0,1,length.out=nrow(tmpSparse$smoothCov)), bw_userCov = tmpSparse$bwCov, rotationCut = c(0, 1), kernel = p2$kernel, rcov = rcovSparse)
sqrt(pccovE_Sparse$sigma2)
