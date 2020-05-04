# devtools::load_all()
#options(error=recover)
library(testthat)

# set up
y <- list(c(1, 2), 4, c(0, 2, 3))
t <- list(c(1.5, 2.5), 2, c(1, 1.5, 2.5))
obsGrid <- seq(0, 3, length.out=7)
mu <- obsGrid
pts <- seq(0, 1, length.out=7)
phi <- cbind(sin(2 * pi * pts), cos(2 * pi * pts))
lambda <- c(6, 1)
rho <- 0.5
fittedCov <- phi %*% diag(lambda) %*% t(phi)


test_that('Truncation works for GetRho', {
  set.seed(1)
  n <- 20
  pts <- signif(seq(0, 1, by=0.05), 14)
  truncPts <- signif(seq(0.1, 0.9, 0.05), 14)
  mu <- rep(0, length(pts))
  samp4 <- Wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
  samp4 <- Sparsify(samp4, pts, 10)
  samp4$Ly[[1]] <- samp4$Lt[[1]] <- c(0, 1)
  samp4Trunc <- TruncateObs(samp4$Ly, samp4$Lt, truncPts)  
  pTrunc <- SetOptions(samp4$Ly, samp4$Lt, list(dataType='Sparse', error=TRUE, kernel='epan', verbose=TRUE))
  smc4 <- GetSmoothedCovarSurface(samp4$Ly, samp4$Lt, mu, pts, pts, pTrunc)
  eig4 <- GetEigenAnalysisResults(smc4$smoothCov, pts, pTrunc)
  phiObs <- ConvertSupport(pts, truncPts, phi=eig4$phi)
  CovObs <- ConvertSupport(pts, truncPts, Cov=eig4$fittedCov)
  
  rho4 <- GetRho(samp4Trunc$Ly, samp4Trunc$Lt, pTrunc, mu[1:length(truncPts)],mu[1:length(truncPts)],truncPts, CovObs, eig4$lambda, phiObs, phiObs, truncPts,smc4$sigma2)
  expect_true(rho4 < 0.2)
})
# # Matlab code:
# y{1} = [1, 2]; y{2} = [4]; y{3} = [0, 2, 3]
# t{1} = [1.5, 2.5]; t{2} = [2]; t{3} = [1, 1.5, 2.5]
# ni = cellfun(@length, y);
# mu = linspace(0, 3, 7);
# out1 = linspace(0, 3, 7);
# pts = linspace(0, 1, 7);
# phi = [sin(2 * pi * pts)', cos(2 * pi * pts)'];
# lambda = [6, 1];
# sigma = 0;
# sig1 = 0.4;
# noeig = 2;
# error = 1;
# method = 'CE';
# shrink = 0;
# regular = 0;
# rho = 0;
# [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular);
# LAMBDA = diag(lambda);
# rho = 0.5;
# subID = [1 3];
# tjID = [2 1 2 3];
# verbose = false;

# getScores2(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink,     out1, regular, muSub, phiSub, LAMBDA, rho, subID, tjID)

# T = range(out1);
# gamma = ((trapz(out1, mu.^2)+sum(lambda))/T)^(0.5);
# alpha = linspace(0.01, 0.22,50);
# rho = gamma*alpha;
# cv_rho(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1,   regular, rho, ni, tjID, verbose)

