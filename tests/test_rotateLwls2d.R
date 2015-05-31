set.seed(1)
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(100, pts, sparsify=2:7)
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=TRUE)

# rotate fit is worse
fit3 <- lwls2d(0.5, kern='epan', xin=rcov3$tpairn, yin=rcov3$cxxn, xout=cbind(pts, pts))
plot(fit3)
sum((fit3 - pts)^2)

val3 <- rotateLwls2d(0.5, 'epan', rcov3$tpairn, rcov3$cxxn, xout=cbind(pts, pts))
plot(val3)
sum((val3 - pts)^2)


# rotate fit is better
gridPoints <- seq(-1, 1, by=0.01)
xin <- expand.grid(gridPoints, gridPoints)
yin <- apply(xin, 1, function(x) 5 * abs(x[1]) + x[2])  + rnorm(nrow(xin), 0, 0.1)
contour(gridPoints, gridPoints, matrix(yin, length(gridPoints)))
datin <- as.data.frame(cbind(rotate(as.matrix(xin), -pi / 4), yin))
names(datin) <- c('x1', 'x2', 'y')
fit1 <- locfit(y ~ lp(x1, x2, h = 1, deg = 2), data = datin, kern = 'epan', maxk=1000)
plot(fit1)
fit2 <- locfit(y ~ lp(x1, x2, h = 1, deg = 1), data = datin, kern = 'epan', maxk=500)
plot(fit2)

plot(yin[xin[, 1] == 0])
val <- predict(fit1, cbind(gridPoints, gridPoints))
plot(val)
val1 <- rotateLwls2d(0.12, 'epan', datin[, 1:2], yin, xout=cbind(gridPoints, gridPoints) / sqrt(2))
plot(val1)


debug(rotateLwls2d)
undebug(rotateLwls2d)
