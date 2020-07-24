### demeanFuc: see FPCReg

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

dx <- function(p, intLenX, gridNumX, brkX){    
	for(i in 1:p){
		if (i == 1) {dxMatrix <- diag(intLenX[1] / (gridNumX[1] - 1), gridNumX[1])}else{dxMatrix <- cdiag(dxMatrix, diag(intLenX[i] / (gridNumX[i]-1), gridNumX[i]))}
		dxMatrix[brkX[i]+1, brkX[i]+1] <- dxMatrix[brkX[i]+1, brkX[i]+1]/2
		dxMatrix[brkX[i+1], brkX[i+1]] <- dxMatrix[brkX[i+1], brkX[i+1]]/2
		}
	dxMatrix <- sqrt(as.matrix(dxMatrix))
	return(dxMatrix)
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

