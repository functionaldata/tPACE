devtools::load_all()
set.seed(1)
n <- 100
pts <- seq(0, 1, by=0.05)
outPts <- seq(0, 1, by=0.1)
samp3 <- wiener(n, pts) + rnorm(n * length(pts), sd=0.5)
samp3 <- sparsify(samp3, pts, 5:10)
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=TRUE)
brcov3 <- BinRawCov(rcov3)

gcv3 <- gcvlwls2d(pts, kern='epan', rcov=brcov3)
# rotate fit is worse
fit3 <- lwls2d(gcv3$h, kern='epan', xin=brcov3$tPairs, yin=brcov3$meanVals, win=brcov3$count, xout=cbind(outPts, outPts))
sum((fit3 - outPts)^2)

val3 <- rotateLwls2d( gcv3$h, 'epan', xin=brcov3$tPairs, yin=brcov3$meanVals, win=brcov3$count, xout=cbind(outPts, outPts))
sum((val3 - outPts)^2)


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
