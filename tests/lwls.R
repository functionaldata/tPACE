library(testthat)
set.seed(1)
x <- seq(0, pi, by=0.1)
y <- sin(x) + rnorm(length(x), sd=0.2)
plot(x, y)
# To MATLAB
write.table(format(data.frame(x=x, y=y), digits=16), file='lwls1.csv', sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)

tmp <- lwls1d(0.3, 'epan', xin=x, yin=y)
plot(tmp)
tmp <- lwls1d(0.3, 'epan', xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3.14))
predict(tmp)

tmp1 <- locfit(y ~ lp(x, h=0.1), data.frame(x, y))
plot(tmp1)
