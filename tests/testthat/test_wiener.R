# source('.../wiener.R')

pts <- seq(0, 1, length=20)
tmp <- wiener(50, pts)
tmp1 <- sparsify(tmp, pts, c(2, 4, 6))

sapply(tmp1$tList, length)
sapply(tmp1$yList, length)

set.seed(1)
tmp2 <- wiener(10, pts=seq(0, 1, by=0.1))
set.seed(1)
tmp3 <- wiener(10, pts=seq(0, 1, by=0.1), sparsify=2)

pts <- seq(0, 1, by=0.02)
tmp <- wiener(1000, pts)
tmp1 <- sparsify(tmp, pts, 1:5, fragment=0.2)
createDesignPlot(tmp1[['tList']], pts, TRUE, FALSE)