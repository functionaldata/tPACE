# devtools::load_all()

viewPdf <- FALSE
set.seed(1)
n <- 100
M <- 51
pts <- seq(0, 1, length.out=M)
mu <- rep(0, length(pts))
sampDense <- Wiener(n, pts)
samp <- Sparsify(sampDense, pts, M)
res <- FPCA(samp$Ly, samp$Lt, list(error=TRUE, FVEthreshold=1, dataType='Dense', plot=TRUE))

if (viewPdf) {
    pdf('tmp.pdf')
}

plot(res)

test_that("inner product of eigenfunctions and mean shall be non-negative", {  
    inprod.mu.and.phi <- res$mu %*% res$phi / M #compute inner product by simple Riemann sum
    expect_true( all(inprod.mu.and.phi>=0) ) #absolute difference
    
})
