## The binned version is exactly the same as the unbinned version.
set.seed(1)
pts <- seq(0, 1, by=0.01)
samp3 <- wiener(2000, pts, sparsify=1:length(pts))
y <- unlist(samp3$yList)
x <- unlist(samp3$tList)
system.time(tmp <- locfit(y ~ lp(x, h=0.1, deg=1)))
plot(tmp)

# bin x:
bins <- sort(unique(x))
system.time({
yBinned <- tapply(y, x, function(yy) c(mean(yy), length(yy)))
yBinned1 <- do.call(rbind, yBinned)
})

system.time(yBinned2 <- aggregate(y, list(x), function(yy) c(mean(yy), length(yy))))

max(abs(as.matrix(yBinned1) - yBinned2[[2]]))

tmp2 <- locfit(yBinned1[, 1] ~ lp(bins, h=0.1, deg=1), weights=yBinned1[, 2])
plot(tmp)
plot(tmp2)
max(abs(predict(tmp, pts) - predict(tmp2, pts)))

