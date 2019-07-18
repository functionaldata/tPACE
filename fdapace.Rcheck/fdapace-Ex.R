pkgname <- "fdapace"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fdapace')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BwNN")
### * BwNN

flush(stderr()); flush(stdout())

### Name: BwNN
### Title: Minimum bandwidth based on kNN criterion.
### Aliases: BwNN

### ** Examples

tinyGrid = list(c(1,7), c(2,3),  6,  c(2,4), c(4,5))
BwNN(tinyGrid, k = 2) # c(3,2)



cleanEx()
nameEx("CreateBWPlot")
### * CreateBWPlot

flush(stderr()); flush(stdout())

### Name: CreateBWPlot
### Title: Functional Principal Component Analysis Bandwidth Diagnostics
###   plot
### Aliases: CreateBWPlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res1 <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=FALSE))
CreateBWPlot(res1)



cleanEx()
nameEx("CreateBasis")
### * CreateBasis

flush(stderr()); flush(stdout())

### Name: CreateBasis
### Title: Create an orthogonal basis of K functions in [0, 1], with nGrid
###   points.
### Aliases: CreateBasis

### ** Examples

basis <- CreateBasis(3, type='fourier')
head(basis)




cleanEx()
nameEx("CreateCovPlot")
### * CreateCovPlot

flush(stderr()); flush(stdout())

### Name: CreateCovPlot
### Title: Create the covariance surface plot based on the results from
###   FPCA() or FPCder().
### Aliases: CreateCovPlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
CreateCovPlot(res)



cleanEx()
nameEx("CreateDesignPlot")
### * CreateDesignPlot

flush(stderr()); flush(stdout())

### Name: CreateDesignPlot
### Title: Create the design plot of the functional data.
### Aliases: CreateDesignPlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
CreateDesignPlot(sampWiener$Lt, sort(unique(unlist(sampWiener$Lt))))



cleanEx()
nameEx("CreateFuncBoxPlot")
### * CreateFuncBoxPlot

flush(stderr()); flush(stdout())

### Name: CreateFuncBoxPlot
### Title: Create functional boxplot using 'bagplot', 'KDE' or 'pointwise'
###   methodology
### Aliases: CreateFuncBoxPlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
CreateFuncBoxPlot(res, list(addIndx=c(1:3)) )



cleanEx()
nameEx("CreateModeOfVarPlot")
### * CreateModeOfVarPlot

flush(stderr()); flush(stdout())

### Name: CreateModeOfVarPlot
### Title: Functional Principal Component Analysis mode of variation plot
### Aliases: CreateModeOfVarPlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
CreateModeOfVarPlot(res)



cleanEx()
nameEx("CreateOutliersPlot")
### * CreateOutliersPlot

flush(stderr()); flush(stdout())

### Name: CreateOutliersPlot
### Title: Functional Principal Component or Functional Singular Value
###   Decomposition Scores Plot using 'bagplot' or 'KDE' methodology
### Aliases: CreateOutliersPlot

### ** Examples

## Not run: 
##D set.seed(1)
##D n <- 420
##D pts <- seq(0, 1, by=0.05)
##D sampWiener <- Wiener(n, pts)
##D sampWiener <- Sparsify(sampWiener, pts, 10)
##D res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
##D             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
##D CreateOutliersPlot(res)
## End(Not run)



cleanEx()
nameEx("CreatePathPlot")
### * CreatePathPlot

flush(stderr()); flush(stdout())

### Name: CreatePathPlot
### Title: Create the fitted sample path plot based on the results from
###   FPCA().
### Aliases: CreatePathPlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan',
            verbose=TRUE))
CreatePathPlot(res, subset=1:5)

# CreatePathPlot has a lot of usages:
## Not run: 
##D CreatePathPlot(res)
##D CreatePathPlot(res, 1:20)
##D CreatePathPlot(res, 1:20, showObs=FALSE)
##D CreatePathPlot(res, 1:20, showMean=TRUE, showObs=FALSE)
##D CreatePathPlot(res, 1:20, obsOnly=TRUE)
##D CreatePathPlot(res, 1:20, obsOnly=TRUE, showObs=FALSE)
##D CreatePathPlot(inputData=sampWiener, subset=1:20, obsOnly=TRUE)
## End(Not run)




cleanEx()
nameEx("CreateScreePlot")
### * CreateScreePlot

flush(stderr()); flush(stdout())

### Name: CreateScreePlot
### Title: Create the scree plot for the fitted eigenvalues
### Aliases: CreateScreePlot

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
CreateScreePlot(res)



cleanEx()
nameEx("CreateStringingPlot")
### * CreateStringingPlot

flush(stderr()); flush(stdout())

### Name: CreateStringingPlot
### Title: Create plots for observed and stringed high dimensional data
### Aliases: CreateStringingPlot

### ** Examples

set.seed(1)
n <- 50
wiener = Wiener(n = n)[,-1]
p = ncol(wiener)
rdmorder = sample(size = p, x=1:p, replace = FALSE)
stringingfit = Stringing(X = wiener[,rdmorder], disOptns = "correlation")
diff_norev = sum(abs(rdmorder[stringingfit$StringingOrder] - 1:p))
diff_rev = sum(abs(rdmorder[stringingfit$StringingOrder] - p:1))
if(diff_rev <= diff_norev){
  stringingfit$StringingOrder = rev(stringingfit$StringingOrder)
  stringingfit$Ly = lapply(stringingfit$Ly, rev)
}
CreateStringingPlot(stringingfit, 1:20)




cleanEx()
nameEx("DynCorr")
### * DynCorr

flush(stderr()); flush(stdout())

### Name: DynCorr
### Title: Dynamical Correlation
### Aliases: DynCorr

### ** Examples

set.seed(10)
n=200             # sample size
t=seq(0,1,length.out=100)       # length of data
mu_quad_x=8*t^2-4*t+5
mu_quad_y=8*t^2-12*t+6
fun=rbind(rep(1,length(t)),-t,t^2)
z1=matrix(0,n,3)
z1[,1]=rnorm(n,0,2)
z1[,2]=rnorm(n,0,16/3)
z1[,3]=rnorm(n,0,4)
x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))
for (i in 1:n){
  x1_quad_error[i,]=mu_quad_x+z1[i,]%*%fun+rnorm(length(t),0,0.01)
  y1_quad_error[i,]=mu_quad_y+2*z1[i,]%*%fun +rnorm(length(t),0,0.01)
}
dyn1_quad=DynCorr(x1_quad_error,y1_quad_error,t) 



cleanEx()
nameEx("Dyn_test")
### * Dyn_test

flush(stderr()); flush(stdout())

### Name: Dyn_test
### Title: Bootstrap test of Dynamic correlation
### Aliases: Dyn_test

### ** Examples

n=200             # sample size
t=seq(0,1,length.out=100)       # length of data
mu_quad_x=8*t^2-4*t+5
mu_quad_y=8*t^2-12*t+6
fun=rbind(rep(1,length(t)),-t,t^2)
z1=matrix(0,n,3)
z1[,1]=rnorm(n,0,2)
z1[,2]=rnorm(n,0,16/3)
z1[,3]=rnorm(n,0,4)   # covariance matrix of random effects
x1_quad_error=y1_quad_error=matrix(0,nrow=n,ncol=length(t))
for (i in 1:n){
  x1_quad_error[i,]=mu_quad_x+z1[i,]%*%fun+rnorm(length(t),0,0.01)
  y1_quad_error[i,]=mu_quad_y+2*z1[i,]%*%fun +rnorm(length(t),0,0.01)
}
bt_DC=Dyn_test(x1_quad_error,y1_quad_error,t,B=1000)




cleanEx()
nameEx("FAM")
### * FAM

flush(stderr()); flush(stdout())

### Name: FAM
### Title: Functional Additive Models
### Aliases: FAM

### ** Examples

set.seed(1000)

library(MASS)

f1 <- function(t) 0.5*t
f2 <- function(t) 2*cos(2*pi*t/4)
f3 <- function(t) 1.5*sin(2*pi*t/4)
f4 <- function(t) 2*atan(2*pi*t/4)

n<-250
N<-500

sig <- diag(c(4.0,2.0,1.5,1.2))

scoreX <- mvrnorm(n,mu=rep(0,4),Sigma=sig)
scoreXTest <- mvrnorm(N,mu=rep(0,4),Sigma=sig)

Y <- f1(scoreX[,1]) + f2(scoreX[,2]) + f3(scoreX[,3]) + f4(scoreX[,4]) + rnorm(n,0,0.1)
YTest <- f1(scoreXTest[,1]) + f2(scoreXTest[,2]) + 
  f3(scoreXTest[,3]) + f4(scoreXTest[,4]) + rnorm(N,0,0.1)

phi1 <- function(t) sqrt(2)*sin(2*pi*t)
phi2 <- function(t) sqrt(2)*sin(4*pi*t)
phi3 <- function(t) sqrt(2)*cos(2*pi*t)
phi4 <- function(t) sqrt(2)*cos(4*pi*t)

grid <- seq(0,1,length.out=21)
Lt <- Lx <- list()
for (i in 1:n) {
  Lt[[i]] <- grid
  Lx[[i]] <- scoreX[i,1]*phi1(grid) + scoreX[i,2]*phi2(grid) + 
    scoreX[i,3]*phi3(grid) + scoreX[i,4]*phi4(grid) + rnorm(1,0,0.01)
}

LtTest <- LxTest <- list()
for (i in 1:N) {
  LtTest[[i]] <- grid
  LxTest[[i]] <- scoreXTest[i,1]*phi1(grid) + scoreXTest[i,2]*phi2(grid) + 
    scoreXTest[i,3]*phi3(grid) + scoreXTest[i,4]*phi4(grid) + rnorm(1,0,0.01)
}


# estimation
fit <- FAM(Y=Y,Lx=Lx,Lt=Lt)

xi <- fit$xi

par(mfrow=c(2,2))
j <- 1
g1 <- f1(sort(xi[,j]))
tmpSgn <- sign(sum(g1*fit$fam[,j]))
plot(sort(xi[,j]),g1,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi1')
points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')

j <- 2
g2 <- f2(sort(xi[,j]))
tmpSgn <- sign(sum(g2*fit$fam[,j]))
plot(sort(xi[,j]),g2,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi2')
points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')

j <- 3
g3 <- f3(sort(xi[,j]))
tmpSgn <- sign(sum(g3*fit$fam[,j]))
plot(sort(xi[,j]),g3,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi3')
points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')

j <- 4
g4 <- f4(sort(xi[,j]))
tmpSgn <- sign(sum(g4*fit$fam[,j]))
plot(sort(xi[,j]),g4,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi4')
points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')


# fitting
fit <- FAM(Y=Y,Lx=Lx,Lt=Lt,nEval=0)
yHat <- fit$mu+apply(fit$fam,1,'sum')
par(mfrow=c(1,1))
plot(yHat,Y)
abline(coef=c(0,1),col=2)


# R^2
R2 <- 1-sum((Y-yHat)^2)/sum((Y-mean(Y))^2)
R2


# prediction
fit <- FAM(Y=Y,Lx=Lx,Lt=Lt,newLx=LxTest,newLt=LtTest)
yHat <- fit$mu+apply(fit$fam,1,'sum')
par(mfrow=c(1,1))
plot(yHat,YTest,xlim=c(-10,10))
abline(coef=c(0,1),col=2)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("FCCor")
### * FCCor

flush(stderr()); flush(stdout())

### Name: FCCor
### Title: Calculate functional correlation between two simultaneously
###   observed processes.
### Aliases: FCCor

### ** Examples

set.seed(1)
n <- 200
nGridIn <- 50
sparsity <- 1:5 # must have length > 1
bw <- 0.2
kern <- 'epan'
T <- matrix(seq(0.5, 1, length.out=nGridIn))

## Corr(X(t), Y(t)) = 1/2
A <- Wiener(n, T)
B <- Wiener(n, T) 
C <- Wiener(n, T) + matrix((1:nGridIn) , n, nGridIn, byrow=TRUE)
X <- A + B
Y <- A + C
indEach <- lapply(1:n, function(x) sort(sample(nGridIn, sample(sparsity, 1))))
tAll <- lapply(1:n, function(i) T[indEach[[i]]])
Xsp <- lapply(1:n, function(i) X[i, indEach[[i]]])
Ysp <- lapply(1:n, function(i) Y[i, indEach[[i]]])

plot(T, FCCor(Xsp, Ysp, tAll, bw)[['corr']], ylim=c(-1, 1))
abline(h=0.5)



cleanEx()
nameEx("FCReg")
### * FCReg

flush(stderr()); flush(stdout())

### Name: FCReg
### Title: Functional Concurrent Regression by 2D smoothing method.
### Aliases: FCReg

### ** Examples

# Y(t) = \beta_0(t) + \beta_1(t) X_1(t) + \beta_2(t) Z_2 + \epsilon

# Settings
set.seed(1)
n <- 75
nGridIn <- 150
sparsity <- 5:10 # Sparse data sparsity
T <- round(seq(0, 1, length.out=nGridIn), 4) # Functional data support
bw <- 0.1
outGrid <- round(seq(min(T), 1, by=0.05), 2)

# Simulate functional data 
mu <- T * 2 # mean function for X_1
sigma <- 1

beta_0 <- 0
beta_1 <- 1
beta_2 <- 1

Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(mu, n, nGridIn, byrow=TRUE)
epsilon <- rnorm(n, sd=sigma)
Y <- matrix(NA, n, nGridIn)
for (i in seq_len(n)) {
  Y[i, ] <- beta_0 + beta_1 * X_1[i, ] + beta_2 * Z[i, 2] + epsilon[i]
}

# Sparsify functional data
set.seed(1)
X_1sp <- Sparsify(X_1, T, sparsity)
set.seed(1)
Ysp <- Sparsify(Y, T, sparsity)
vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
withError2D <- FCReg(vars, bw, bw, outGrid)



cleanEx()
nameEx("FClust")
### * FClust

flush(stderr()); flush(stdout())

### Name: FClust
### Title: Functional clustering and identifying substructures of
###   longitudinal data
### Aliases: FClust

### ** Examples

## Not run: 
##D data(medfly25)
##D Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs) 
##D newClust <- FClust(Flies$Ly, Flies$Lt, k = 2, optnsFPCA = 
##D                     list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90))
##D                     
##D # We denote as 'veryLowCount' the group of flies that lay less
##D # than twenty-five eggs during the 25-day period examined.
##D 
##D veryLowCount = ifelse( sapply( unique(medfly25$ID), function(u) 
##D                    sum( medfly25$nEggs[medfly25$ID == u] )) < 25, 0, 1)
##D N <- length(unique(medfly25$ID))
##D (correctRate <- sum( (1 + veryLowCount) ==  newClust$cluster) / N) # 99.6%
## End(Not run)



cleanEx()
nameEx("FOptDes")
### * FOptDes

flush(stderr()); flush(stdout())

### Name: FOptDes
### Title: Optimal Designs for Functional and Longitudinal Data for
###   Trajectory Recovery or Scalar Response Prediction
### Aliases: FOptDes

### ** Examples

set.seed(1)
n <- 50
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- MakeFPCAInputs(IDs = rep(1:n, each=length(pts)), 
                             tVec = rep(pts, times = n), 
                             yVec = t(sampWiener))
res <- FOptDes(Ly=sampWiener$Ly, Lt=sampWiener$Lt, p=2,
               isSequential=FALSE, RidgeCand = seq(0.02,0.2,0.02))



cleanEx()
nameEx("FPCA")
### * FPCA

flush(stderr()); flush(stdout())

### Name: FPCA
### Title: Functional Principal Component Analysis
### Aliases: FPCA

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
plot(res) # The design plot covers [0, 1] * [0, 1] well.
CreateCovPlot(res, 'Fitted')



cleanEx()
nameEx("FPCAder")
### * FPCAder

flush(stderr()); flush(stdout())

### Name: FPCAder
### Title: Take derivative of an FPCA object
### Aliases: FPCAder

### ** Examples


bw <- 0.2
kern <- 'epan'
set.seed(1)
n <- 100 
M <- 40
pts <- seq(0, 1, length.out=M)
lambdaTrue <- c(1, 0.8, 0.1)^2
sigma2 <- 0.1

samp2 <- MakeGPFunctionalData(n, M, pts, K=length(lambdaTrue), 
                              lambda=lambdaTrue, sigma=sqrt(sigma2), basisType='legendre01')
samp2 <- c(samp2, MakeFPCAInputs(tVec=pts, yVec=samp2$Yn))
fpcaObj <- FPCA(samp2$Ly, samp2$Lt, list(methodMuCovEst='smooth',
                userBwCov=bw, userBwMu=bw, kernel=kern, error=TRUE)) 
CreatePathPlot(fpcaObj, showObs=FALSE)

FPCoptn <- list(bw=bw, kernelType=kern, method='FPC')
DPCoptn <- list(bw=bw, kernelType=kern, method='DPC')
FPC <- FPCAder(fpcaObj, FPCoptn)
DPC <- FPCAder(fpcaObj, DPCoptn)

CreatePathPlot(FPC, ylim=c(-5, 10))
CreatePathPlot(DPC, ylim=c(-5, 10))

# Get the true derivatives
phi <-  CreateBasis(K=3, type='legendre01', pts=pts)
basisDerMat <- apply(phi, 2, function(x) 
                       ConvertSupport(seq(0, 1, length.out=M - 1), pts, diff(x) * (M - 1)))
trueDer <- matrix(1, n, M, byrow=TRUE) + tcrossprod(samp2$xi, basisDerMat)
matplot(t(trueDer), type='l', ylim=c(-5, 10))

# DPC is slightly better in terms of RMSE
mean((fitted(FPC) - trueDer)^2)
mean((fitted(DPC) - trueDer)^2)




cleanEx()
nameEx("FVPA")
### * FVPA

flush(stderr()); flush(stdout())

### Name: FVPA
### Title: Functional Variance Process Analysis for dense functional data
### Aliases: FVPA

### ** Examples

set.seed(1)
n <- 25
pts <- seq(0, 1, by=0.01)
sampWiener <- Wiener(n, pts)
# Data have to dense for FVPA to be relevant!
sampWiener <- Sparsify(sampWiener, pts, 101) 
fvpaObj <- FVPA(sampWiener$Ly, sampWiener$Lt)



cleanEx()
nameEx("GetCrCovYX")
### * GetCrCovYX

flush(stderr()); flush(stdout())

### Name: GetCrCovYX
### Title: Functional Cross Covariance between longitudinal variable Y and
###   longitudinal variable X
### Aliases: GetCrCovYX

### ** Examples

Ly1= list( rep(2.1,7), rep(2.1,3),2.1 );
Lt1 = list(1:7,1:3, 1);
Ly2 = list( rep(1.1,7), rep(1.1,3),1.1); 
Lt2 = list(1:7,1:3, 1);
Ymu1 = rep(55,7);
Ymu2 = rep(1.1,7);
AA<-GetCrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2)
  



cleanEx()
nameEx("GetCrCovYZ")
### * GetCrCovYZ

flush(stderr()); flush(stdout())

### Name: GetCrCovYZ
### Title: Functional Cross Covariance between longitudinal variable Y and
###   scalar variable Z
### Aliases: GetCrCovYZ

### ** Examples

Ly <- list( runif(5),  c(1:3), c(2:4), c(4))
Lt <- list( c(1:5), c(1:3), c(1:3), 4)
Z = rep(4,4) # Constant vector so the covariance has to be zero.
sccObj = GetCrCovYZ(bw=1, Z= Z, Ly=Ly, Lt=Lt, Ymu=rep(4,5))



cleanEx()
nameEx("GetNormalisedSample")
### * GetNormalisedSample

flush(stderr()); flush(stdout())

### Name: GetNormalisedSample
### Title: Normalise sparse functional sample
### Aliases: GetNormalisedSample GetNormalizedSample

### ** Examples

set.seed(1)
n <- 100
M <- 51
pts <- seq(0, 1, length.out=M)
mu <- rep(0, length(pts))
sampDense <- MakeGPFunctionalData(n, M, mu, K=1, basisType='sin', sigma=0.01)
samp4 <- MakeFPCAInputs(tVec=sampDense$pts, yVec=sampDense$Yn)
res4E <- FPCA(samp4$Ly, samp4$Lt, list(error=TRUE))
sampN <- GetNormalisedSample(res4E, errorSigma=TRUE)

CreatePathPlot(subset=1:20, inputData=samp4, obsOnly=TRUE, showObs=FALSE)
CreatePathPlot(subset=1:20, inputData=sampN, obsOnly=TRUE, showObs=FALSE)



cleanEx()
nameEx("MultiFAM")
### * MultiFAM

flush(stderr()); flush(stdout())

### Name: MultiFAM
### Title: Functional Additive Models with Multiple Predictor Processes
### Aliases: MultiFAM

### ** Examples

set.seed(1000)

library(MASS)

f11 <- function(t) t
f12 <- function(t) 2*cos(2*pi*t/4)
f21 <- function(t) 1.5*sin(2*pi*t/4)
f22 <- function(t) 1.5*atan(2*pi*t/4)

n<-100
N<-200

sig <- matrix(c(2.0, 0.0, 0.5, -.2,
                0.0, 1.2, -.2, 0.3,
                0.5, -.2, 1.7, 0.0,
                -.2, 0.3, 0.0, 1.0),
              nrow=4,ncol=4)

scoreX <- mvrnorm(n,mu=rep(0,4),Sigma=sig)
scoreXTest <- mvrnorm(N,mu=rep(0,4),Sigma=sig)

Y <- f11(scoreX[,1]) + f12(scoreX[,2]) + f21(scoreX[,3]) + f22(scoreX[,4]) + rnorm(n,0,0.5)
YTest <- f11(scoreXTest[,1]) + f12(scoreXTest[,2]) + 
f21(scoreXTest[,3]) + f22(scoreXTest[,4]) + rnorm(N,0,0.5)

phi11 <- function(t) sqrt(2)*sin(2*pi*t)
phi12 <- function(t) sqrt(2)*sin(4*pi*t)
phi21 <- function(t) sqrt(2)*cos(2*pi*t)
phi22 <- function(t) sqrt(2)*cos(4*pi*t)

grid <- seq(0,1,length.out=21)
Lt <- Lx1 <- Lx2 <- list()
for (i in 1:n) {
  Lt[[i]] <- grid
  Lx1[[i]] <- scoreX[i,1]*phi11(grid) + scoreX[i,2]*phi12(grid) + rnorm(1,0,0.01)
  Lx2[[i]] <- scoreX[i,3]*phi21(grid) + scoreX[i,4]*phi22(grid) + rnorm(1,0,0.01)
}

LtTest <- Lx1Test <- Lx2Test <- list()
for (i in 1:N) {
  LtTest[[i]] <- grid
  Lx1Test[[i]] <- scoreXTest[i,1]*phi11(grid) + scoreXTest[i,2]*phi12(grid) + rnorm(1,0,0.01)
  Lx2Test[[i]] <- scoreXTest[i,3]*phi21(grid) + scoreXTest[i,4]*phi22(grid) + rnorm(1,0,0.01)
}

X1 <- list(Ly=Lx1, Lt=Lt)
X2 <- list(Ly=Lx2, Lt=Lt)

X1Test <- list(Ly=Lx1Test, Lt=LtTest)
X2Test <- list(Ly=Lx2Test, Lt=LtTest)

X <- list(X1, X2)
XTest <- list(X1Test, X2Test)

# estimation
sbf <- MultiFAM(Y=Y,X=X)

xi <- sbf$xi

par(mfrow=c(2,2))
j <- 1
p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
g11 <- f11(sort(xi[,j])) - 
trapzRcpp(sort(xi[,j]),f11(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
tmpSgn <- sign(sum(g11*sbf$SBFit[,j]))
plot(sort(xi[,j]),g11,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi11')
points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)

j <- 2
p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
g12 <- f12(sort(xi[,j])) - 
trapzRcpp(sort(xi[,j]),f12(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
tmpSgn <- sign(sum(g12*sbf$SBFit[,j]))
plot(sort(xi[,j]),g12,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi12')
points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)

j <- 3
p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
g21 <- f21(sort(xi[,j])) - 
trapzRcpp(sort(xi[,j]),f21(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
tmpSgn <- sign(sum(g21*sbf$SBFit[,j]))
plot(sort(xi[,j]),g21,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi21')
points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)

j <- 4
p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
g22 <- f22(sort(xi[,j])) - 
trapzRcpp(sort(xi[,j]),f22(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
tmpSgn <- sign(sum(g22*sbf$SBFit[,j]))
plot(sort(xi[,j]),g22,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi22')
points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)

## Not run: 
##D # fitting
##D sbf <- MultiFAM(Y=Y,X=X,nEval=0)
##D yHat <- sbf$mu+apply(sbf$SBFit,1,'sum')
##D par(mfrow=c(1,1))
##D plot(yHat,Y)
##D abline(coef=c(0,1),col=2)
##D 
##D 
##D # R^2
##D R2 <- 1-sum((Y-yHat)^2)/sum((Y-mean(Y))^2)
##D R2
##D 
##D 
##D # prediction
##D sbf <- MultiFAM(Y=Y,X=X,XTest=XTest)
##D yHat <- sbf$mu+apply(sbf$SBFit,1,'sum')
##D par(mfrow=c(1,1))
##D plot(yHat,YTest)
##D abline(coef=c(0,1),col=2)
## End(Not run)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("SBFitting")
### * SBFitting

flush(stderr()); flush(stdout())

### Name: SBFitting
### Title: Iterative Smooth Backfitting Algorithm
### Aliases: SBFitting

### ** Examples

set.seed(100)

n <- 100
d <- 2
X <- pnorm(matrix(rnorm(n*d),nrow=n,ncol=d)%*%matrix(c(1,0.6,0.6,1),nrow=2,ncol=2))

f1 <- function(t) 2*(t-0.5)
f2 <- function(t) sin(2*pi*t)

Y <- f1(X[,1])+f2(X[,2])+rnorm(n,0,0.1)

# component function estimation
N <- 101
x <- matrix(rep(seq(0,1,length.out=N),d),nrow=N,ncol=d)
h <- c(0.12,0.08)
  
sbfEst <- SBFitting(Y,x,X,h)
fFit <- sbfEst$SBFit

par(mfrow=c(1,2))
plot(x[,1],f1(x[,1]),type='l',lwd=2,col=2,lty=4,xlab='X1',ylab='Y')
points(x[,1],fFit[,1],type='l',lwd=2,col=1)
points(X[,1],Y,cex=0.3,col=8)
legend('topleft',legend=c('SBF','true'),col=c(1,2),lwd=2,lty=c(1,4),horiz=FALSE,bty='n')
abline(h=0,col=8)

plot(x[,2],f2(x[,2]),type='l',lwd=2,col=2,lty=4,xlab='X2',ylab='Y')
points(x[,2],fFit[,2],type='l',lwd=2,col=1)
points(X[,2],Y,cex=0.3,col=8)
legend('topright',legend=c('SBF','true'),col=c(1,2),lwd=2,lty=c(1,4),horiz=FALSE,bty='n')
abline(h=0,col=8)

# prediction
x <- X
h <- c(0.12,0.08)
  
sbfPred <- SBFitting(Y,X,X,h)
fPred <- sbfPred$mY+apply(sbfPred$SBFit,1,'sum')

par(mfrow=c(1,1))
plot(fPred,Y,cex=0.5,xlab='SBFitted values',ylab='Y')
abline(coef=c(0,1),col=8)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("Stringing")
### * Stringing

flush(stderr()); flush(stdout())

### Name: Stringing
### Title: Stringing for High-Dimensional data
### Aliases: Stringing

### ** Examples

set.seed(1)
n <- 50
wiener = Wiener(n = n)[,-1]
p = ncol(wiener)
rdmorder = sample(size = p, x=1:p, replace = FALSE)
stringingfit = Stringing(X = wiener[,rdmorder], disOptns = "correlation")
diff_norev = sum(abs(rdmorder[stringingfit$StringingOrder] - 1:p))
diff_rev = sum(abs(rdmorder[stringingfit$StringedOrder] - p:1))
if(diff_rev <= diff_norev){
  stringingfit$StringingOrder = rev(stringingfit$StringingOrder)
  stringingfit$Ly = lapply(stringingfit$Ly, rev)
}
plot(1:p, rdmorder[stringingfit$StringingOrder], pch=18); abline(a=0,b=1)




cleanEx()
nameEx("WFDA")
### * WFDA

flush(stderr()); flush(stdout())

### Name: WFDA
### Title: Warped Functional DAta Analysis
### Aliases: WFDA

### ** Examples

N = 44;
eps = 0.123;
M = 41;
set.seed(123) 
Tfinal = 3
me <- function(t) exp(-Tfinal*(((t/Tfinal^2)-0.5))^2);
T = seq(0,Tfinal,length.out = M) 
recondingTimesMat = matrix(nrow = N, ncol = M)
yMat = matrix(nrow = N, ncol = M)

for (i in 1:N){
  peak = runif(min = 0.2,max =  0.8,1) * Tfinal 
  recondingTimesMat[i,] = sort( unique(c( seq(0.0 , peak, length.out = round((M+1)/2)),
                            seq( peak, Tfinal, length.out = round((M+1)/2))))) * Tfinal
  yMat[i,] = me(recondingTimesMat[i,]) * rnorm(1, mean=4.0, sd=  eps)
                                       + rnorm(M, mean=0.0, sd=  eps) 
}

Y = as.list(as.data.frame(t(yMat)))
X = rep(list(T),N)
 
sss =  WFDA(Ly = Y, Lt = X, list( choice = 'weighted' ))
par(mfrow=c(1,2))
matplot(x= T, t(yMat), t='l', main = 'Raw', ylab = 'Y'); grid()
matplot(x= T, t(sss$aligned), t='l', main = 'Aligned', ylab='Y'); grid() 



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fitted.FPCA")
### * fitted.FPCA

flush(stderr()); flush(stdout())

### Name: fitted.FPCA
### Title: Fitted functional sample from FPCA object
### Aliases: fitted.FPCA

### ** Examples

set.seed(1)
n <- 100
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 5:10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
fittedY <- fitted(res, ciOptns = list(alpha=0.05))

workGrid <- res$workGrid
cvgUpper <- fittedY$cvgUpper
cvgLower <- fittedY$cvgLower

par(mfrow=c(2,3))
ind <- sample(1:n,6)
for (i in 1:6) {
 j <- ind[i]
 plot(workGrid,cvgUpper[j,],type='l',ylim=c(min(cvgLower[j,]),max(cvgUpper[j,])),col=4,lty=2,
   xlab='t', ylab='X(t)', main=paste(j,'-th subject',sep=''))
 points(workGrid,cvgLower[j,],type='l',col=4,lty=2)
 points(res$inputData$Lt[[j]],res$inputData$Ly[[j]])
}
    



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fitted.FPCAder")
### * fitted.FPCAder

flush(stderr()); flush(stdout())

### Name: fitted.FPCAder
### Title: Fitted functional sample from FPCAder object
### Aliases: fitted.FPCAder

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)



cleanEx()
nameEx("kCFC")
### * kCFC

flush(stderr()); flush(stdout())

### Name: kCFC
### Title: Functional clustering and identifying substructures of
###   longitudinal data using kCFC.
### Aliases: kCFC

### ** Examples

data(medfly25) 
Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs)
kcfcObj <- kCFC(Flies$Ly[1:250], Flies$Lt[1:250], # using only 250 for speed consideration 
                 optnsSW = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90),
                 optnsCS = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.70))



cleanEx()
nameEx("plot.FPCA")
### * plot.FPCA

flush(stderr()); flush(stdout())

### Name: CreateDiagnosticsPlot
### Title: Functional Principal Component Analysis Diagnostics plot
### Aliases: CreateDiagnosticsPlot plot.FPCA

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res1 <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', error=FALSE, kernel='epan', verbose=FALSE))
plot(res1)



cleanEx()
nameEx("predict.FPCA")
### * predict.FPCA

flush(stderr()); flush(stdout())

### Name: predict.FPCA
### Title: Predict FPC scores for a new sample given an FPCA object
### Aliases: predict.FPCA

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt)
res




cleanEx()
nameEx("print.FPCA")
### * print.FPCA

flush(stderr()); flush(stdout())

### Name: print.FPCA
### Title: Print an FPCA object
### Aliases: print.FPCA

### ** Examples

set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- Sparsify(sampWiener, pts, 10)
res <- FPCA(sampWiener$Ly, sampWiener$Lt)
res




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
