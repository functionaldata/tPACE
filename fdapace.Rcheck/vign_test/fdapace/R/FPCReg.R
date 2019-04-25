#' Function for performing functonal linear regression where the covariates are functions X1(t1),X2(t2),.. and the response is a function Y(t_y).
#'
#' @param vars A list of input functional covariates with name of "X1", "X2",.. and a functional response with name "Y". Each field should have two fields: 'Lt', a list (sparse) or a matrix (Dense) specifying the time of observations, and 'Ly', a list (Sparse) or a matrix (Dense) of the observations.
#' @param varsOptns A list of options named by "X1", "X2",..."Y". Each filed specify the paramaters that control the corresponding variables. (default: see details of FPCA())
#' @param isNewSub A 1*n vector of 0s or 1s, where n is the total count of subjects. 0 denotes the corresponding subject is only used for estimation and 1 denotes the corresponding subject is only used for prediction. (default: 0's)
#' @param method The method used for selecting the number of principal components of functional predictors X's used in functional regression , including 'AIC', 'BIC' and 'FVE'. (default: "AIC")
#' @param FVEthreshold A scalar specifying the proportion used for 'FVE'. (default: 0.99)
#' @param alpha A scalar specifying the level of the confidence bands. (default: 0.05)
#' @param Kx The number of principal components of functional predictors X's used in functional regression.
#'
#' @return A list containing the following fields:
#' \item{estiBeta}{A list with fields of estimated beta_XiY(s,t) defiend on [range(Xi),range(Y)]}
#' \item{predictY}{A list containing fitted or predicted (when is NewSub is true) functions for E(Y|X).}
#' \item{cbandY}{A list with confidence bands of E(Y|X).}
#' \item{Q}{Quasi R-square}
#' \item{r2}{Functional R-square.}
#' \item{varsMean}{A list with mean function of covariates and response.}
#' \item{Kx}{The number of principal components of functional predictors X's used in functional regression.}
#'
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.}
#'
#' @export
#' @examples
#' set.seed(1000)
#' #Model: E(Y(t)|X) = int(beta(s,t)*X(s))
#' n <- 200 #number of subjects
#' ngrids <- 51 #number of grids in [0,1] for X(s)
#' ngridt <- 101 #number of grids in [0,1] for Y(t)
#' grids <- seq(0, 1, length.out=ngrids) #regular grids in [0,1] for X(s)
#' gridt <- seq(0, 1, length.out=ngridt) #regular grids in [0,1] for Y(t)
#' 
#' #generate X
#' #{1, sqrt(2)*sin(2*pi*s), sqrt(2)*cos(2*pi*t)} are used to generate X.
#' eigenFun <- list( function(s){1 + 0 * s}, 
#'                   function(s){sqrt(2) * sin(2*pi*s)},
#'                   function(s){sqrt(2) * cos(2*pi*s)})
#' 
#' sig <- matrix(c(1.5, 0.0, 0.0, 0.9, -.5, 0.1,
#'                 0.0, 1.2, 0.0, -.3, 0.8, 0.4,
#'                 0.0, 0.0, 1.0, 0.4, -.3, 0.7,
#'                 0.9, -.3, 0.4, 2.0, 0.0, 0.0,
#'                 -.5, 0.8, -.3, 0.0, 1.5, 0.0,
#'                 0.1, 0.4, 0.7, 0.0, 0.0, 1.0),
#'                 nrow=6,ncol=6)
#' 
#' scoreX <- MASS::mvrnorm(n,mu=rep(0,6),Sigma=sig)
#' scoreX1 <- scoreX[,1:3]
#' scoreX2 <- scoreX[,4:6]
#' 
#' basisX1 <- sapply(eigenFun,function(x){x(grids)})
#' latentX1 <- scoreX1 %*% t(basisX1)
#' measErrX1 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
#' denseX1 <- latentX1 + measErrX1
#' 
#' basisX2 <- sapply(eigenFun,function(x){x(grids)})
#' latentX2 <- scoreX2 %*% t(basisX2)
#' measErrX2 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
#' denseX2 <- latentX2 + measErrX2
#' 
#' #generate Y
#' #beta(s, t) <- sin(2 * pi * s)*cos(2 * pi * t)
#' betaEigen1 <- function(t){f <- function(s){
#'                             sin(2*pi*s) * cos(2*pi*t) * (1+0*s)}; return(f)}
#' betaEigen2 <- function(t){f <- function(s){
#'                             sin(2*pi*s) * cos(2*pi*t) * (sqrt(2)*sin(2*pi*s))}; return(f)}
#' betaEigen3 <- function(t){f <- function(s){
#'                             sin(2*pi*s) * cos(2*pi*t) * (sqrt(2)*cos(2*pi*s))}; return(f)}
#' betaEigen <- list(betaEigen1, betaEigen2, betaEigen3) 
#' basisY <- array(0,c(ngridt, 3))
#' for(i in 1:3){
#' 	intbetaEigen <- function (t) {integrate(betaEigen[[i]](t), lower = 0, upper = 1)$value}
#' 	basisY[, i] <- sapply(1:ngridt, function(x){intbetaEigen(gridt[x])})
#' 	}
#' latentY <- scoreX1 %*% t(basisY) - scoreX2 %*% t(basisY)
#' measErrY <- sqrt(0.01) * matrix(rnorm(n*ngridt), n, ngridt) #0.01 is sigma^2
#' denseY <- latentY + measErrY
#' 
#' #======Dense data===============================================
#' timeX <- t(matrix(rep(grids, n),length(grids), n))
#' timeY <- t(matrix(rep(gridt, n),length(gridt), n))
#' denseVars <- list(X1 = list(Ly = denseX1, Lt = timeX),
#'                   X2 = list(Ly = denseX2, Lt = timeX),
#'                   Y=list(Ly = denseY,Lt = timeY))
#' 
#' resuDense <- FPCReg(denseVars, method="FVE") 
#' 
#' par(mfrow=c(1,2))
#' estiBetaX1Y_Dense <- resuDense$estiBeta$betaX1Y
#' args1 <- list(xlab = 's', ylab = 't', zlab = 'estiBetaX1Y_Dense(s, t)',
#'               lighting = FALSE, phi = 45, theta = 45)
#' args2 <- list(x = 1:ngrids, y = 1:ngridt, z = estiBetaX1Y_Dense[1:ngrids, 1:ngridt])
#' do.call(plot3D::persp3D,c(args2, args1))
#' 
#' estiBetaX2Y_Dense <- resuDense$estiBeta$betaX2Y
#' args1 <- list(xlab = 's', ylab = 't', zlab = 'estiBetaX2Y_Dense(s, t)',
#'              lighting = FALSE, phi = 45, theta = 45)
#' args2 <- list(x = 1:ngrids, y = 1:ngridt, z = estiBetaX2Y_Dense[1:ngrids, 1:ngridt])
#'  # do.call(plot3D::persp3D,c(args2, args1))
#' 
#' #======Sparse data===============================================
#' \dontrun{
#' sparsity = 5:8
#' sparseX1 <- Sparsify(denseX1, grids, sparsity)
#' sparseX2 <- Sparsify(denseX2, grids, sparsity)
#' sparseY <- Sparsify(denseY, gridt, sparsity)
#' sparseVars <- list(X1 = sparseX1, X2 = sparseX2, Y = sparseY)
#' 
#' resuSparse <- FPCReg(sparseVars, method="FVE", FVEthreshold=0.98) 
#' #or resuSparse <- FPCReg(vars = sparseVars,
#' #                        varsOptns = list(X1=list(userBwCov = 0.03)))
#' 
#' par(mfrow=c(1,2))
#' estiBetaX1Y_Sparse = resuSparse$estiBeta$betaX1Y
#' args1 = list(xlab = 's', ylab = 't', zlab = 'estiBetaX1Y_Sparse(s,t)', 
#'              lighting = FALSE, phi = 45,theta = 45)
#' args2 = list(x = 1:51, y = 1:51, z = estiBetaX1Y_Sparse[1:51, 1:51])
#' do.call(plot3D::persp3D, c(args2, args1))
#' 
#' estiBetaX2Y_Sparse = resuSparse$estiBeta$betaX2Y
#' args1 = list(xlab = 's', ylab = 't', zlab = 'estiBetaX2Y_Sparse(s,t)', 
#'              lighting = FALSE, phi = 45,theta = 45)
#' args2 = list(x = 1:51, y = 1:51, z = estiBetaX2Y_Sparse[1:51, 1:51])
#' do.call(plot3D::persp3D, c(args2, args1))
#' 
#' par(mfrow=c(2,3))
#' for(i in 1:6){
#' 	plot(sparseVars[['Y']]$Lt[[i]], sparseVars[['Y']]$Ly[[i]], 
#' 	xlab = 'time', ylab = 'observations', ylim = c(-1.5, 1.5))
#' 	lines(seq(0, 1, length.out = 51), resuSparse$predictY[[i]])
#' 	lines(seq(0, 1, length.out = 51), resuSparse$cbandY[[i]][,2], lty = 2)
#' 	lines(seq(0, 1, length.out = 51), resuSparse$cbandY[[i]][,1], lty = 2)
#' 	}
#' 	}

FPCReg <- function(vars, varsOptns = NULL, isNewSub = NULL, method = 'AIC', FVEthreshold = 0.99, alpha = 0.05, Kx = NULL){
	#===============data checking and manipulation===============
	p <- length(vars) - 1
	if (p == 0) stop('Too few covariates.')
	if (!'Y' %in% names(vars)|"" %in% names(vars)) stop('Missing name of the response which should be "Y".')
	if (!'X1' %in% names(vars)|"" %in% names(vars)) stop('Missing name of preictors which should be "X1","X2",...')
	if (anyDuplicated(names(vars)) > 0) stop('Duplicated Names.')
	if ('Y' %in% names(vars)) {vars <- c(vars[names(vars) != 'Y'], vars['Y'])}
	
	#Set Lt in [0,1] for dense data without Lt. 
	for (i in 1:(p+1)) {
		if (is.null(vars[[i]]$Lt) & is.matrix(vars[[i]]$Ly)) {
			n1 <- dim(vars[[i]]$Ly)[1]
			n2 <- dim(vars[[i]]$Ly)[2]
			vars[[i]]$Lt <- t(matrix(rep(seq(0, 1, length.out = n2), n1), n2, n1))
		}
	}
	
	#Dense==1 iff all covariates and response are Matrix
	if (sum(sapply(vars, function(x){is.matrix(x$Ly)}) * sapply(vars, function(x){is.matrix(x$Lt)}))) {Dense <- 1} else {Dense <- 0} 
	
	#set options 
	#varsOptnsDef is the default
	#specific default setting for optns now and output grids for sparse is fixed at 51 now due to the cross-cov.
	if (Dense == 1) {optns <- list(dataType = "Dense" ,error = 1, kernel='gauss', nRegGrid=51, useBinnedData='OFF')} else {optns <- list(dataType = "Sparse", error = TRUE, kernel = 'gauss' ,nRegGrid = 51, useBinnedData = 'OFF')} 
	for (i in 1:(p+1)) {
		if (i == 1) {varsOptnsDef <- list(optns)} else {varsOptnsDef <- c(varsOptnsDef, list(optns))}
	}
	names(varsOptnsDef) <- names(vars)
	if (is.null(varsOptns)) {varsOptns <- varsOptnsDef} else if (!is.list(varsOptns)) {stop('wrong inpiut of varsOptns.')}
	if (!sum(names(varsOptns) %in% names(vars)) == length(names(varsOptns))) stop('Check names of varsOptns which should be X1,X2..Y.')
	for (i in 1:(p+1)) {
		name1=names(varsOptns)
		if (!names(vars)[i]%in%names(varsOptns)) {varsOptns <- c(varsOptns, list(optns)); names(varsOptns) = c(name1, names(vars)[i])}
		if (is.null(varsOptns[names(vars)[i]]$dataType) & Dense == 0) {varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], dataType = "Sparse")
		}else if(is.null(varsOptns[names(vars)[i]]$dataType) & Dense==1){
		varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], dataType = "Dense")
		}
		if (is.null(varsOptns[names(vars)[i]]$nRegGrid) & Dense == 0) {varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], nRegGrid = 51)}
		if (is.null(varsOptns[names(vars)[i]]$error)) {varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], error = TRUE)}
		if (Dense == 1) {varsOptns[[names(vars)[i]]]$nRegGrid <- dim(vars[[i]]$Ly)[2]}
	}
	varsOptns <- varsOptns[names(vars)]
	
	#nRegGrids determine the numbers of ouput grids 
	if (Dense == 0) {nRegGrids <- sapply(varsOptns, function(x){x$nRegGrid})}

	for ( i in 1:(p+1)) {
		if (!'Ly' %in% names(vars[[i]])) stop('Insert the name "Ly" for the predictors and response to indicate time of data.')
		if (!'Lt' %in% names(vars[[i]])) stop('Insert the name "Lt" for the predictors and response to indicate time of data.')
      	if (! length(names(vars[[i]])) == 2) stop('Check data.')
	}

	#list dense data for using HandleNumericsAndNAN and demeanFuc func.
	for (i in 1:(p+1)) {
		if (is.list(vars[[i]]$Lt) == 0) {vars[[i]]$Lt <- lapply(1:nrow(vars[[i]]$Lt), function(j) vars[[i]]$Lt[j, ])}
		if (is.list(vars[[i]]$Ly) == 0) {vars[[i]]$Ly <- lapply(1:nrow(vars[[i]]$Ly), function(j) vars[[i]]$Ly[j, ])}
	}
	vars[sapply(vars, is.list)] <- lapply(vars[sapply(vars, is.list)], function(v) HandleNumericsAndNAN(v[['Ly']], v[['Lt']]))

	if (is.null(isNewSub)) {isNewSub <- rep(0, length(vars[[1]]$Lt))}
	varsTrain <- vars
	for (i in 1:(p+1)) {
		varsTrain[[i]]$Lt <- vars[[i]]$Lt[which(isNewSub == 0)]
		varsTrain[[i]]$Ly <- vars[[i]]$Ly[which(isNewSub == 0)]	
	}

	#===============population parameters===============
	demeanedRes <- demeanFuc(p, varsTrain, kern='gauss', varsOptns) #Centered predictors. Using gauss for demeanFuc, but may be relaxed.
	varsTrain <- demeanedRes[['xList']]
	muList <- demeanedRes[['muList']]
	intLen <- array(0, p+1) #lengthes of time window 
	gridNum <- array(0, p+1) #number of grids for varss
	regGrid <- list()
	varsMean <- list()
	for (i in 1:(p+1)) {
		intLen[i] <- max(c(unlist(varsTrain[[i]]$Lt)))-min(c(unlist(varsTrain[[i]]$Lt)))	
		gridNum[i] <- if (Dense==1) {length(varsTrain[[i]]$Ly[[1]])} else {nRegGrids[i]}      #sparse output grids is fixed to be 51 so far.
		regGrid[[i]] <- seq(min(unlist(varsTrain[[i]]$Lt)), max(unlist(varsTrain[[i]]$Lt)), length.out=gridNum[i]) #time grids for output mean and beta. 
		varsMean[[i]] <- muList[[i]](regGrid[[i]]) 
	}
	brkX <- c(0, cumsum(gridNum[1:p]))      #break points for matrix products
	dxMatrix <- dx(p, intLen[1:p], gridNum[1:p], brkX) #generate the inner product for multivariate predictors
	
	#matrilize the dense data.
	if (Dense == 1) {
		for (i in 1:(p+1)) {
			varsTrain[[i]]$Lt <- do.call(rbind, varsTrain[[i]]$Lt)
			varsTrain[[i]]$Ly <- do.call(rbind, varsTrain[[i]]$Ly)
		}
	}

	varsCov <- if (Dense == 1) {denseCov(p, varsTrain, brkX, dxMatrix, gridNum, varsOptns)} else {sparseCov(p, varsTrain, brkX, varsOptns, vars, muList, gridNum, dxMatrix, regGrid)}
	xCov <- varsCov[['xCov']]
	croCov <- varsCov[['croCov']]
	yCov <- varsCov[['yCov']]
	varsSigma2 <- varsCov[['varsSigma2']]
	numposiEigen <- varsCov[['numposiEigen']] #number of positive eigenvalues
	eigenValue <- varsCov[['eigenValue']]	
	eigenFun <- varsCov[['eigenFun']]
	
	#estimate score
	if (Dense == 1) {
		diagSigma2 <- varsCov[['diagSigma2']]
		for (i in 1:p) {
			if (i == 1) {varsMatrix <- t(sapply(vars[[1]]$Ly, c))}else{varsMatrix <- cbind(varsMatrix, t(sapply(vars[[i]]$Ly, c)))}
			if (i == 1) {muListX <- muList[[1]](regGrid[[i]])}else{muListX <- c(muListX, muList[[i]](regGrid[[i]]))}
		}
		varsMatrixDemean <- apply(varsMatrix,1, '-', muListX)
		varsMatrixDemean <- t(varsMatrixDemean)
		denseScore <- innerProd(varsMatrixDemean, eigenFun, brkX, intLen[1:p]) #for dense data, he score is directly estimated by integral(x,ei)
	} else {
		varsDemean <- varsCov[['varsDemean']]
		subCov <- varsCov[['subCov']]	
		subEigenFun <- varsCov[['subEigenFun']]	
		subMeanX <- varsCov[['subMeanX']]	
		diagSigma2 <- varsCov[['diagSigma2']]	
		numposiEigen <- varsCov[['numposiEigen']]
	    #for sparse data, empute score is estimated in sparseCov()
		invSigma <- lapply(subCov, function(x){as.matrix(solve(x))})
		invSigmaY <- mapply(function(X, Y){X %*% Y}, X=invSigma, Y=varsDemean, SIMPLIFY=FALSE)
		emputeScore <- mapply(function(X, Y){diag(eigenValue[1:numposiEigen]) %*% X %*% Y}, X = subEigenFun, Y = invSigmaY, SIMPLIFY = FALSE)
	}

	#if (Dense == 0) {Kx <- sparseComp(emputeScore, varsDemean, method, FVEthreshold, eigenValue, numposiEigen, isNewSub, diagSigma2, subEigenFun)} else {Kx <- denseComp(method, varsMatrixDemean, denseScore, eigenFun, eigenValue, isNewSub, diagSigma2, numposiEigen, FVEthreshold)}
	if(!is.null(Kx)==1) {Kx <- Kx} else if (Dense == 0) {Kx <- sparseComp(emputeScore, varsDemean, method, FVEthreshold, eigenValue, numposiEigen, isNewSub, diagSigma2, subEigenFun)} else {Kx <- denseComp(method, varsMatrixDemean, denseScore, eigenFun, eigenValue, isNewSub, diagSigma2, numposiEigen, FVEthreshold)}

	
	#prodedure of FPCA Regression
	pcaBeta <- array(0, c(dim(croCov)))
	for (i in 1:Kx) {
		pcaBeta <- pcaBeta + (1/eigenValue[i]) * eigenFun[, i] %*% innerProd(t(eigenFun[, i]), croCov, brkX, intLen)
	}
	estiBeta <- list()
	#divide beta for output
	for (i in 1:p) {
		estiBeta[[i]] <- pcaBeta[(brkX[i]+1):brkX[i+1],]
	}
	names(estiBeta) <- lapply(names(vars)[1:p], function(x){paste("beta", x, "Y", sep = '')})

	intBetaPhi <- innerProd(t(pcaBeta), as.matrix(eigenFun[, 1:Kx]), brkX, intLen)
	diagYcov <- diag(yCov) ; diagYcov[which(diagYcov<0)] <- min(diagYcov[diagYcov>0])
	varY <- innerProd(t(diagYcov), as.matrix(rep(1, gridNum[p+1])), c(0,gridNum[p+1]), intLen[p+1])
	varYX <- sum(innerProd(t(intBetaPhi^2), as.matrix(rep(1, gridNum[p+1])),c(0, gridNum[p+1]), intLen[p+1]) * eigenValue[1:Kx])

	Q <- min(varYX / varY, 1) #quasi R^2?
    
	r2 <- apply(intBetaPhi^2 %*% diag(eigenValue[1:Kx]), 1, sum) / diagYcov #functional R^2? 
	r2 <- sapply(r2, function(x){min(x, 1)})
	#R2=innerProd(t(r2),as.matrix(rep(1,gridNum[p+1])),c(0,gridNum[p+1]),intLen[p+1])/intLen[p+1]

	#===============the predictions===============
	if (Dense == 1) {
		meanY <- varsMean[[p+1]]
		names(varsMean) <- names(vars)
		predictY <- meanY + intBetaPhi%*%t(denseScore[, 1:Kx])
		predictY <- t(predictY)
		diageValue <- diag(eigenValue[1:Kx])
		for (i in 1:p) {
			if (i == 1) {
				diagSigma2 <- c(rep(varsSigma2[1], gridNum[1]))
			}else{
				diagSigma2 = c(diagSigma2, rep(varsSigma2[i],gridNum[i]))
			}
		}
		SigmaY <- xCov + diag(diagSigma2)  #This part may be considered to use raw cov
		Omega <- diageValue - diageValue %*% t(eigenFun[, 1:Kx]) %*% solve(SigmaY) %*% eigenFun[, 1:Kx] %*% t(diageValue)
		pwVar <- diag(intBetaPhi %*% Omega %*% t(intBetaPhi))
		upCI <- t(apply(predictY, 1, '+', stats::qnorm(1-alpha / 2) * sqrt(pwVar)))
		lwCI <- t(apply(predictY, 1, '-', stats::qnorm(1-alpha / 2) * sqrt(pwVar)))
		cbandY <- list(upCI = upCI, lwCI = lwCI)
	} else {
		score <- lapply(emputeScore, function(x){x[1:Kx]})		
		imputeY <- lapply(score, function(x){intBetaPhi %*% x})
		meanY <- varsMean[[p+1]]
		names(varsMean) <- names(vars)
		predmean <- innerProd(t(pcaBeta), as.matrix(subMeanX), brkX, intLen[1:p])
		predictY <- lapply(imputeY, function(x){x + meanY})
		diageValue <- diag(eigenValue[1:Kx])
		H <- lapply(subEigenFun, function(x){diag(eigenValue[1:Kx]) %*% x[1:Kx, ]})
		Omega <- mapply(function(X, Y){diageValue - X %*% Y %*% t(X)}, X=H, Y=invSigma, SIMPLIFY = FALSE)
		pwVar <- lapply(Omega, function(x){diag(intBetaPhi %*% x %*% t(intBetaPhi))})
		upCI <- mapply(function(X, Y){X + stats::qnorm(1 - alpha / 2)*sqrt(Y)}, X = predictY, Y = pwVar)
		lwCI <- mapply(function(X, Y){X - stats::qnorm(1 - alpha / 2)*sqrt(Y)}, X = predictY, Y = pwVar)
		cbandY <- lapply(1:dim(lwCI)[2], function(x){CI = cbind(lwCI[, x], upCI[, x]) ; colnames(CI) = c('lwCI', 'upCI') ; return(CI)})
	}

	res <- list(estiBeta = estiBeta, predictY = predictY, cbandY = cbandY, Q = Q, r2 = r2, Kx = Kx, varsMean = varsMean)
	res
}

demeanFuc <- function(p, varsTrain, kern, varsOptns) {
	for (i in 1:(p+1)) {
	    userBwMuXi <- SetOptions(varsTrain[[i]]$Ly, varsTrain[[i]]$Lt, varsOptns[[i]])$userBwMu
		varsTrain[[i]] <- list(Lt = varsTrain[[i]]$Lt, Ly = varsTrain[[i]]$Ly, userBwMu = userBwMuXi)
	}
	tmp <- lapply(varsTrain, function(x) {
		Tin <- sort(unique(unlist(x[['Lt']])))
		xmu <- GetSmoothedMeanCurve(x[['Ly']], x[['Lt']], Tin, Tin[1], list(userBwMu = x[['userBwMu']], kernel=kern))[['mu']]
		muFun<-approxfun(Tin, xmu)
		x[['Ly']] <- lapply(1:length(x[['Ly']]), function(i) x[['Ly']][[i]] - muFun(x[['Lt']][[i]]))
		xmu <- muFun
		list(x = x, mu = xmu)
  	})
  	xList <- lapply(tmp, `[[`, 'x')
  	muList <- lapply(tmp, `[[`, 'mu')
  	list(xList = xList, muList = muList)
}

innerProd <- function(A, B, gridbrk, intlen){      
	brkX <- gridbrk
	oupmat <- array(0, c(dim(A)[1], dim(B)[2]))
	for (i in 1:(length(brkX)-1)) {
		a1 <- brkX[i] + 1
		a2 <- brkX[i + 1]
		A[,a1] <- A[, a1] / 2 ; A[, a2] = A[, a2] / 2
		oupmat <- oupmat + A[, a1:a2] %*% B[a1:a2, ] * intlen[i] / (a2 - a1)
	}
	return(oupmat)
}

dx <- function(p, intLenX, gridNumX, brkX){    
	for(i in 1:p){
		if (i == 1) {dxMatrix <- diag(intLenX[1] / (gridNumX[1] - 1), gridNumX[1])}else{dxMatrix <- cdiag(dxMatrix, diag(intLenX[i] / (gridNumX[i]-1), gridNumX[i]))}
		dxMatrix[brkX[i]+1, brkX[i]+1] <- dxMatrix[brkX[i]+1, brkX[i]+1]/2
		dxMatrix[brkX[i+1], brkX[i+1]] <- dxMatrix[brkX[i+1], brkX[i+1]]/2
	}
	dxMatrix <- sqrt(as.matrix(dxMatrix))
	return(dxMatrix)
}

denseCov <- function(p, varsTrain, brkX, dxMatrix, gridNum, varsOptns){
	indRaw <- varsTrain[[1]]$Ly
	if (varsOptns[[1]]$error == 1) {covDense <- GetCovDense(indRaw, optns=varsOptns[[1]]) ; blkCov <- covDense$smoothCov ; varsSigma2 = covDense$sigma2} else {blkCov <- cov(indRaw)}
	if (p >= 2) {
		for (i in 2:p) {	
			indRaw <- varsTrain[[i]]$Ly		
			if (varsOptns[[i]]$error == 1) {covDense <- GetCovDense(indRaw, optns=varsOptns[[i]]) ; blkCov <- cdiag(blkCov, covDense$smoothCov) ; varsSigma2 <- c(varsSigma2, covDense$sigma2)} else {blkCov <- cdiag(blkCov, cov(indRaw))}	
		}

		for (i in 1:p) {
			gridBrk1 <- brkX[i]+1
			gridBrk2 <- brkX[i+1]
			for (j in 1:p) {
				if (j > i) {
					gridBrk3 <- brkX[j]+1
					gridBrk4 <- brkX[j+1]
					blkCov[gridBrk1:gridBrk2, gridBrk3:gridBrk4] <- cov(varsTrain[[i]]$Ly, varsTrain[[j]]$Ly)
					blkCov[gridBrk3:gridBrk4, gridBrk1:gridBrk2] <- t(as.matrix(blkCov[gridBrk1:gridBrk2, gridBrk3:gridBrk4]))
				}
			}
		}
	}
	xCov <- as.matrix(blkCov)
	adjxCov <- eigen(dxMatrix %*% xCov %*% dxMatrix)
	eigenValue <- adjxCov$value
	eigenValue <- sapply(eigenValue, function(x){max(x, 0)})
	numposiEigen <- sum(eigenValue > 0)
	eigenFun <- diag(diag(1 / dxMatrix)) %*% adjxCov$vector
	xCov <- eigenFun %*% diag(eigenValue) %*% t(eigenFun)

	for (i in 1:p) {
		if(i==1) {croCov <- cov(varsTrain[[i]]$Ly, varsTrain[[p+1]]$Ly)} else {croCov <- rbind(croCov, cov(varsTrain[[i]]$Ly, varsTrain[[p+1]]$Ly))}
	}
	indRaw <- varsTrain[[p+1]]$Ly
	if (varsOptns[[p+1]]$error == 1) {yCov <- GetCovDense(indRaw, optns = varsOptns[[p+1]])$smoothCov} else {yCov <- cov(indRaw)}
	for (i in 1:p) {
		if(i == 1){diagSigma2 <- c(rep(varsSigma2[1], gridNum[1]))}else{diagSigma2 <- c(diagSigma2, rep(varsSigma2[i], gridNum[i]))}
	}
	return(list(xCov = xCov, croCov = croCov, yCov = yCov, varsSigma2 = varsSigma2, eigenValue = eigenValue, eigenFun = eigenFun, diagSigma2 = diagSigma2, numposiEigen = numposiEigen))
}

#gridNum may be used in the future to generalize the output grids, not restricted by 51.
sparseCov <- function(p, varsTrain, brkX, varsOptns, vars, muList, gridNum, dxMatrix, regGrid){
	Ly <- varsTrain[[1]]$Ly ; Lt=varsTrain[[1]]$Lt;
	optns <- SetOptions(Ly, Lt, varsOptns[[1]])
	obsGrid <- sort(unique(c(unlist(Lt))))
	meanCr <- array(0,c(length(obsGrid)))
	meanCrAll <- list(meanCr)    #meanCrAll is a input required by GetCrCovYX(). It's zero since varsTrain is demeaned.
	rfpca <- FPCA(Ly, Lt, optns)
	blkCov <- rfpca$smoothedCov
	varsSigma2 <- rfpca$sigma2
	varsbwMu <- list(rfpca$bwMu)
	varsbwCov <- list(rfpca$bwCov)

	if (p>=2) {
		for (i in 2:p) {
			Ly <- varsTrain[[i]]$Ly ; Lt = varsTrain[[i]]$Lt;
			optns <- SetOptions(Ly, Lt, varsOptns[[i]])
			obsGrid <- sort(unique(c(unlist(Lt))))
			meanCr <- array(0, c(length(obsGrid)))
			meanCrAll[[i]] <- meanCr
			rfpca <- FPCA(Ly, Lt, optns)
			mx_temp1 <- rfpca$smoothedCov
			assign(paste("sigma_", i, sep = ""), rfpca$sigma2)
			varsSigma2 <- c(varsSigma2, get(paste("sigma_", i, sep = "")))
			blkCov <- cdiag(blkCov, mx_temp1)	
			varsbwMu[[i]] <- rfpca$bwMu
			varsbwCov[[i]] <- rfpca$bwCov
		}

		for (i in 1:p) {
			gridBrk1 <- brkX[i]+1
			gridBrk2 <- brkX[i+1]
			for (j in 1:p) {
				if (j>i) {
					gridBrk3 <- brkX[j]+1
					gridBrk4 <- brkX[j+1]
					blkCov[gridBrk1:gridBrk2, gridBrk3:gridBrk4] <- GetCrCovYX(bw1 = varsbwCov[[i]], bw2 = varsbwCov[[j]], Ly1 = varsTrain[[i]]$Ly, Lt1 = varsTrain[[i]]$Lt, Ymu1 = meanCrAll[[i]], Ly2 = varsTrain[[j]]$Ly, Lt2 = varsTrain[[j]]$Lt, Ymu2 = meanCrAll[[j]])$smoothedCC
					blkCov[gridBrk3:gridBrk4, gridBrk1:gridBrk2] <- t(as.matrix(blkCov[gridBrk1:gridBrk2, gridBrk3:gridBrk4]))
				}
			}
		}
	}
	xCov <- as.matrix(blkCov)
	adjxCov <- eigen(dxMatrix %*% xCov %*% dxMatrix)
	eigenValue <- adjxCov$value
	eigenValue <- sapply(eigenValue, function(x){max(x, 0)})
	numposiEigen <- sum(eigenValue > 0)
	eigenFun <- diag(diag(1 / dxMatrix)) %*% adjxCov$vector
	xCov <- eigenFun %*% diag(eigenValue) %*% t(eigenFun)
	
	Ly <- varsTrain[[p+1]]$Ly ; Lt <- varsTrain[[p+1]]$Lt;
	optns <- SetOptions(Ly, Lt, varsOptns[[p+1]])
	obsGrid <- sort(unique(c(unlist(Lt))))
	yfpca <- FPCA(Ly, Lt, optns)
	yCov <- yfpca$fittedCov
	varsbwMu[[p+1]] <- yfpca$bwMu
	varsbwCov[[p+1]] <- yfpca$bwCov
	meanCrY <- array(0, c(length(obsGrid)))
	for(i in 1:p){
		if (i==1) {croCov <- GetCrCovYX(bw1 = varsbwCov[[i]], bw2 = varsbwCov[[p+1]], Ly1 = varsTrain[[i]]$Ly, Lt1 = varsTrain[[i]]$Lt, Ymu1 = meanCrAll[[i]], Ly2 = varsTrain[[p+1]]$Ly, Lt2 = varsTrain[[p+1]]$Lt, Ymu2 = meanCrY)$smoothedCC
		}else {
		croCov <- rbind(croCov,GetCrCovYX(bw1 = varsbwCov[[i]], bw2 = varsbwCov[[p+1]], Ly1 = varsTrain[[i]]$Ly, Lt1 = varsTrain[[i]]$Lt, Ymu1 = meanCrAll[[i]], Ly2 = varsTrain[[p+1]]$Ly, Lt2 = varsTrain[[p+1]]$Lt, Ymu2 = meanCrY)$smoothedCC)
		}
	}
	#==========================
	for (i in 1:p) {
		varDemean <- mapply('-',vars[[i]]$Ly, lapply(vars[[i]]$Lt, muList[[i]]), SIMPLIFY = FALSE)
		if(i==1){varsDemean <- varDemean}else{varsDemean <- mapply(function(X, Y){c(X, Y)}, X = varsDemean, Y = varDemean, SIMPLIFY = FALSE)}
	}
	nZero <- rep(list(0), lengthVarsFPCReg(vars))
	for (i in 1:p) {
		varPnumber <- sapply(vars[[i]]$Lt, length)
		nZero <- mapply(c, nZero, varPnumber, SIMPLIFY=FALSE)
	}
	varsPnumber <- lapply(nZero, cumsum)
	if(p == 1){diagSigma2 <- lapply(nZero, function(x){if(x[2]==1){sv <- as.matrix(rep(varsSigma2[1], x[2]))}else{sv <- rep(varsSigma2[1], x[2])};C <- diag(sv)})
	}else{
		diagSigma2 <- lapply(nZero, function(x){if(x[2]==1){sv <- as.matrix(rep(varsSigma2[1], x[2]))}else{sv <- rep(varsSigma2[1], x[2])}; C <- diag(sv) ; for (i in 2:p) {if(x[i+1]==1){sv <- as.matrix(rep(varsSigma2[i], x[i+1]))}else{sv <- rep(varsSigma2[i], x[i+1])};C <- cdiag(C, diag(sv))} ; return(as.matrix(C))})
	}
	
	a1 <- brkX[1]+1
	a2 <- brkX[1+1]
	subCov <- lapply(vars[[1]]$Lt,function(x){ConvertSupport(regGrid[[1]],x,Cov= xCov[a1:a2, a1:a2])})

	if (p>=2) {
		for (i in 2:p) {
			a1 <- brkX[i]+1
			a2 <- brkX[i+1]
			assign(paste("subMx_", i, i, sep = ""),lapply(vars[[i]]$Lt,function(x){ConvertSupport(regGrid[[i]],x,Cov= xCov[a1:a2, a1:a2])}))	
			subCov <- mapply(function(X, Y){cdiag(X, Y)}, subCov, get(paste("subMx_", i, i, sep = "")),SIMPLIFY=FALSE)	
			}
	
		for (i in 1:p) {
			gridBrk1 <- brkX[i]+1
			gridBrk2 <- brkX[i+1]
			for(j in 1:p){
				if (j > i) {
					gridBrk3 <- brkX[j]+1
					gridBrk4 <- brkX[j+1]
					assign(paste("subMx_", i, j, sep = ""), mapply(function(X, Y){gd <- expand.grid(X,Y); matrix(interp2lin(regGrid[[i]], regGrid[[j]], xCov[gridBrk1:gridBrk2, gridBrk3:gridBrk4], gd$Var1, gd$Var2), nrow=length(X))}, X = vars[[i]]$Lt, Y = vars[[j]]$Lt,SIMPLIFY=FALSE))
					subMxx <- get(paste("subMx_", i, j,sep = ""))
					subCov <- mapply(function(X, Y, Z){
						pnumber1 = Z[i] + 1
						pnumber2 = Z[i + 1]
						pnumber3 = Z[j] + 1
						pnumber4 = Z[j + 1]
						X[pnumber1:pnumber2, pnumber3:pnumber4] = Y
						X[pnumber3:pnumber4, pnumber1:pnumber2] = t(Y)
						return(X)
					}, X = subCov, Y = subMxx, Z = varsPnumber ,SIMPLIFY=FALSE)
				}
			}
		}
	}
	subCov <- mapply('+', subCov, diagSigma2, SIMPLIFY=FALSE)
	listEigenFun <- lapply(1:numposiEigen, function(i) eigenFun[,i])

	for (i in 1:p) {
		a1 <- brkX[i]+1
		a2 <- brkX[i+1]
		eigenFunXi <- eigenFun[a1:a2]
		listEigenFunXi <- lapply(listEigenFun,'[', a1:a2)
		phimx0 <- lapply(vars[[i]]$Lt,function(Y){t(sapply(listEigenFunXi, function(x){approxfun(regGrid[[i]], x)(Y)}))})
		phimx0 <- lapply(phimx0,function(x){if(dim(x)[1]==1){x <- t(x)}else{x <- x}})
		if (i == 1) {subEigenFun <- phimx0} else {subEigenFun <- mapply(function(X, Y){cbind(X, Y)}, X = subEigenFun, Y = phimx0 ,SIMPLIFY=FALSE)}
		if (i == 1) {subMeanX <- muList[[i]](regGrid[[i]])} else {subMeanX <- c(subMeanX, muList[[i]](regGrid[[i]]))}
	}

	return(list(xCov = xCov, croCov = croCov, varsSigma2 = varsSigma2, yCov = yCov, varsbwMu = varsbwMu, varsbwCov = varsbwCov, varsDemean = varsDemean, subCov = subCov, subEigenFun = subEigenFun, subMeanX = subMeanX, numposiEigen = numposiEigen, diagSigma2 = diagSigma2, eigenValue = eigenValue, eigenFun = eigenFun))
}

lengthVarsFPCReg <- function(varsTrain, subset) {
	lenEach <- sapply(varsTrain, function(x) {
		if (is.list(x)) {
		x[['userBwMu']] <- NULL
		x[['userBwCov']] <- NULL
			sapply(x, length)
		} else {
			stop('Cannot subset variable')
			}
		}, simplify = FALSE)
  	len <- unique(unlist(lenEach))
  	if (length(len) != 1) {
		stop('Length of variables are not the same!')
	}
	return(len)
}

sparseComp <- function(emputeScore, varsDemean, method, FVEthreshold, eigenValue, numposiEigen, isNewSub, diagSigma2, subEigenFun){
	if (length(which(isNewSub == 1)) == 0) {varsDemeanTrain <- varsDemean}else{varsDemeanTrain <- varsDemean[-which(isNewSub == 1)]}
	if (length(which(isNewSub == 1)) == 0) {emputeScoreTrain <- emputeScore}else{emputeScoreTrain <- emputeScore[-which(isNewSub == 1)]}
	if (length(which(isNewSub == 1)) == 0) {subEigenFunTrain <- subEigenFun}else{subEigenFunTrain <- subEigenFun[-which(isNewSub == 1)]}
	if (length(which(isNewSub == 1)) == 0) {diagSigma2Train <- diagSigma2}else{diagSigma2Train <- diagSigma2[-which(isNewSub == 1)]}
	if (method == 'FVE') {ratio <- cumsum(eigenValue[which(eigenValue > 0)]) / sum(eigenValue[which(eigenValue > 0)]); Kx <- min(which(ratio > FVEthreshold))}	
	if (method == 'AIC') {
		psuLogli <- array(0, c(numposiEigen))
		for (i in 1:numposiEigen) {
			if(i == 1) {psuVec <- mapply(function(X, Y, Z){Z - as.matrix(X[1:i, ]) * Y[1:i]}, X = subEigenFunTrain, Y = emputeScoreTrain, Z = varsDemeanTrain, SIMPLIFY=FALSE)} else {psuVec <- mapply(function(X, Y, Z){Z - as.matrix(t(X[1:i, ]))%*%as.matrix(Y[1:i])}, X = subEigenFunTrain, Y = emputeScoreTrain, Z = varsDemeanTrain, SIMPLIFY=FALSE)}		
			#psuLogli[i] <- sum(mapply(function(X, Y){t(Y) %*% diag(1/diag(X)) %*% Y / 2}, X = diagSigma2Train, Y = psuVec)) + i
			psuLogli[i] <- sum(mapply(function(X, Y){t(Y) %*% solve(X) %*% Y / 2}, X = diagSigma2Train, Y = psuVec)) + i
		}
		diffRatio <- (max(psuLogli) - psuLogli) / diff(range(psuLogli))
		Kx <- min(which(diffRatio > 0.95))
		#diffRatio=abs(psuLogli/min(psuLogli))
		#Kx=min(which(diffRatio<1.2))
		#Kx=which.min(psuLogli>min(psuLogli)*1.1)
	}

	if (method == 'BIC') {
		psuLogli <- array(0, c(numposiEigen))
		for (i in 1:numposiEigen) {
			if (i==1) {psuVec <- mapply(function(X, Y, Z){Z - as.matrix(X[1:i, ]) * Y[1:i]}, X = subEigenFunTrain, Y = emputeScoreTrain, Z = varsDemeanTrain, SIMPLIFY=FALSE)} else {psuVec <- mapply(function(X, Y, Z){Z - as.matrix(t(X[1:i, ])) %*% as.matrix(Y[1:i])}, X = subEigenFunTrain, Y = emputeScoreTrain, Z = varsDemeanTrain, SIMPLIFY=FALSE)}		
			#psuLogli[i] <- sum(mapply(function(X, Y){t(Y) %*% diag(1 / diag(X)) %*% Y / 2}, X = diagSigma2Train,Y=psuVec)) + i * log(length(unlist(varsDemeanTrain))) / 2
			psuLogli[i] <- sum(mapply(function(X, Y){t(Y) %*% solve(X) %*% Y / 2}, X = diagSigma2Train,Y=psuVec)) + i * log(length(unlist(varsDemeanTrain))) / 2
		}
		diffRatio <- (max(psuLogli) - psuLogli)/diff(range(psuLogli))
		Kx <- min(which(diffRatio > 0.95))
		#diffRatio=abs(psuLogli/min(psuLogli))
		#Kx=min(which(diffRatio<1.2))
		#Kx=which.min(psuLogli>min(psuLogli)*1.2)
	}
	Kx
}


denseComp <- function(method, varsMatrixDemean, denseScore, eigenFun, eigenValue, isNewSub, diagSigma2, numposiEigen, FVEthreshold){
	if (method=='FVE') {ratio <- cumsum(eigenValue[which(eigenValue > 0)]) / sum(eigenValue[which(eigenValue > 0)]) ; Kx <- min(which(ratio > FVEthreshold))}
	if (method=='AIC') {
		if (length(which(isNewSub == 1)) == 0) {varsMatrixDemeanTrain <- varsMatrixDemean} else {varsMatrixDemeanTrain <- varsMatrixDemean[-which(isNewSub == 1), ]}
		if (length(which(isNewSub == 1)) == 0) {denseScoreTrain <- denseScore} else {denseScoreTrain <- denseScore[-which(isNewSub == 1), ]}
		psuLogli <- array(0, c(numposiEigen))
		for (i in 1:numposiEigen) {
			psuVec <- varsMatrixDemeanTrain-denseScoreTrain[, 1:i] %*% t(eigenFun[, 1:i])	
			psuLogli[i] <- sum(diag(psuVec %*% diag(1 / diagSigma2) %*% t(psuVec))) / 2 + i
			}
		diffRatio <- (max(psuLogli) - psuLogli) / diff(range(psuLogli))
		Kx <- min(which(diffRatio > 0.95))
		#diffRatio=abs(psuLogli/min(psuLogli))
		#Kx=min(which(diffRatio<1.3))
	}
	if (method == 'BIC') {
		if (length(which(isNewSub == 1)) == 0) {varsMatrixDemeanTrain <- varsMatrixDemean} else {varsMatrixDemeanTrain <- varsMatrixDemean[-which(isNewSub == 1), ]}
		if (length(which(isNewSub == 1)) == 0) {denseScoreTrain <- denseScore} else {denseScoreTrain <- denseScore[-which(isNewSub == 1), ]}
		total_p <- length(varsMatrixDemeanTrain)
		psuLogli <- array(0, c(numposiEigen))
		for (i in 1:numposiEigen) {
			psuVec <- varsMatrixDemeanTrain - denseScoreTrain[, 1:i] %*% t(eigenFun[, 1:i])	
			psuLogli[i] <- sum(diag(psuVec %*% diag(1 / diagSigma2) %*% t(psuVec))) / 2 + i * log(total_p) / 2
			}
		diffRatio <- (max(psuLogli) - psuLogli) / diff(range(psuLogli))
		Kx <- min(which(diffRatio > 0.95))
		#diffRatio=abs(psuLogli/min(psuLogli))
		#Kx=min(which(diffRatio<1.3))
	}
	return(Kx)
}

cdiag <- function(A,B){
	if(is.matrix(A)==0){A <- as.matrix(A)}
	if(is.matrix(B)==0){B <- as.matrix(B)}
	nrow <- dim(A)[1]+dim(B)[1]
	ncol <- dim(A)[2]+dim(B)[2]
	C <- array(0,c(nrow,ncol))
	C[1:dim(A)[1],1:dim(A)[2]] <- A
	C[(dim(A)[1]+1):(dim(A)[1]+dim(B)[1]),(dim(A)[2]+1):(dim(A)[2]+dim(B)[2])] <- B
	return(C)
}



