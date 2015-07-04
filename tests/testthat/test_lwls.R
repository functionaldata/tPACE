library(testthat)
set.seed(1)
x <- seq(0, pi, by=0.1)
y <- sin(x) + rnorm(length(x), sd=0.2)
#plot(x, y)
# To MATLAB
# write.table(format(data.frame(x=x, y=y), digits=16), file='lwls1.csv', sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
## Matlab code
# addpath('~/PACE_matlab/release2.17/PACE')
# dat = csvread('lwls1.csv');
# [~, val] = lwls(0.3, 'epan', 0, 1, 0, dat(:, 1)', dat(:, 2), ones(1, size(dat, 1)), [0, 0.15, 0.35, 0.55, 3], 0)

# tmp <- lwls1d(0.3, 'epan', xin=x, yin=y)
# plot(tmp)


test_that('Epanechnikov kernel matches with different polynomial degrees', {
  expect_equal(Rlwls1d(0.3, 'epan', xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3), win = rep(1, length(x))), 
    c(-0.075034417378891, 0.142526970448573, 0.380805485912617,   0.573996912896636, 0.232449699086771))
  expect_equal(Rlwls1d(0.3, 'epan', npoly=2, xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3), win = rep(1, length(x))), 
    c(-0.125290762148466, 0.130962171169946, 0.442387233004433, 0.496519823276647, 0.346408319249383))

})

test_that('Epanech. kernel matches different bandwidth choice', {
  expect_equal(Rlwls1d(0.5, 'epan', xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3), win = rep(1, length(x))), 
    c(-0.114913056431249,  0.132079727835864,   0.365292140767686, 0.569502194595946,   0.218632400355849))
  expect_equal(Rlwls1d(0.1, 'epan', xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3), win = rep(1, length(x))), 
    c(-0.125290762148466,  0.084052844902148,   0.534948131880310, 0.488730072830666,   0.412855918365678))
})

test_that('Rectangular kernel matches MATLAB output', {
  expect_equal(Rlwls1d(0.3, 'rect', xin=x, yin=y, xout=c(0.15, 0.35, 0.55, 3), win = rep(1, length(x))),
    c(0.140580458121389,   0.369243683204375,   0.611958566485952, 0.220769364459782))
})

test_that('Derivative works matches MATLAB output', {
  expect_equal(Rlwls1d(0.3, 'epan', nder=1, xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3), win = rep(1, length(x))),
    c(0.985197231385895,   1.767923028863129,   0.837595175365137, 0.929221260582535,  -0.689889054107760))
  expect_equal(Rlwls1d(0.3,'epan', npoly=2, nder=1, xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3), win = rep(1, length(x))),
    c(4.452885020486605,  1.715294894713705,   0.837595175365137, 0.929221260582535,  -1.758251118132247))
})

test_that('Differnt weights match MATLAB output', {
  win <- rep(1, length(x))
  win[1:3] <- 2
  expect_equal(Rlwls1d(0.3, 'epan', xin=x, yin=y, win=win, xout=c(0, 0.15, 0.35, 0.55, 3) ),
    c(-0.075034417378891,  0.123584297605571,   0.353323309535711, 0.573996912896636,   0.232449699086771))
})

test_that('Normal kernel matches', {
  expect_equal(Rlwls1d(0.3, 'gauss', xin=x, yin=y, win = rep(1, length(x)), xout=c(0, 0.15, 0.35, 0.55, 3)), 
    c(-0.062129832105935, 0.134704121324048,   0.361071058656017,   0.553453878429684, 0.207951702919279))
})

x <- seq(0, pi,length.out=12)
y <- sin(x) + (x);


test_that('Normal Variant kernel matches', {
  expect_equal(Rlwls1d(0.3,'gausvar', npoly=1, nder=1, xin=x, yin=y, xout=c( min(x), max(x)), win=rep(1, length(x))),
    c(2.005328714335409, -0.005328714335407))
})

test_that('Quartic kernel matches', {
  expect_equal(Rlwls1d(0.3,'quar', npoly=1, nder=1, xin=x, yin=y, xout=c( min(x), max(x)), win=rep(1, length(x))),
    c( 1.986460839127102,   0.013539160872897 ))
})










##  Normal kernel DO NOT MATCH
# test_that('Normal kernel matches', {
    # expect_equal(lwls1d(0.3, 'gauss', xin=x, yin=y, xout=c(0, 0.15, 0.35, 0.55, 3)), 
                 # c(-0.062129832105935, 0.134704121324048,   0.361071058656017,   0.553453878429684, 0.207951702919279))
# })
if (1==2) { #This ain't gonna work any more
test_that('Can return fit object', {
    expect_equal(predict(lwls1d(0.3, 'epan', xin=x, yin=y, returnFit=TRUE)), 
                 predict(locfit(y ~ lp(x, h=0.3, deg=1), data.frame(x=x, y=y), kern='epan', ev=dat())))
    expect_equal(predict(lwls1d(0.3, 'epan', xin=x, yin=y, xout=c(0, 1, 3), returnFit=TRUE)), 
                 predict(locfit(y ~ lp(x, h=0.3, deg=1), data.frame(x=x, y=y), kern='epan', ev=c(0, 1, 3))))
})
}

# test when h is large the implementations are the same.
