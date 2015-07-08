## The binned version is exactly the same as the unbinned version.
set.seed(1)
pts <- seq(0, 1, by=0.01)
samp3 <- wiener(2000, pts, sparsify=1:length(pts))
y <- unlist(samp3$yList)
x <- unlist(samp3$tList)
