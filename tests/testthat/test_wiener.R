# source('.../Wiener.R')

pts <- seq(0, 1, length=20)
tmp <- Wiener(50, pts)
tmp1 <- Sparsify(tmp, pts, c(2, 4, 6))

sapply(tmp1$Lt, length)
sapply(tmp1$Ly, length)

set.seed(1)
tmp2 <- Wiener(10, pts=seq(0, 1, by=0.1))
set.seed(1)
tmp3 <- Wiener(10, pts=seq(0, 1, by=0.1), sparsify=2)

pts <- seq(0, 1, by=0.02)
tmp <- Wiener(1000, pts)
tmp1 <- Sparsify(tmp, pts, 1:5, fragment=0.2)
CreateDesignPlot(tmp1[['Lt']], pts, TRUE, FALSE)