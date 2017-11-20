
AICchoice<-function(varsTrain,CrMatInfo,EigenInfo,BIC = FALSE){
	ncovariate = length(varsTrain)
	nsubj = length(varsTrain[[1]]$Ly)
	neigen = length(EigenInfo$values)
	score = lapply(CrMatInfo$FPCAlist,function(x)x$xiEst)
	Sigma2 = lapply(CrMatInfo$FPCAlist,function(x)x$sigma2)
	MaxL = min(unlist(lapply(score,ncol)))
	Phi = EigenInfo$vectors[,1:MaxL]
	UnitLen = nrow(Phi)/ncovariate
	LoadingGrid = CrMatInfo$FPCAlist[[1]]$workGrid
	##computing the AIC score
	psuLogli = rep(0,MaxL)
	for(L in 1:MaxL){
			for(i in 1:nsubj){
			U = lapply(varsTrain,function(x)x$Ly[[i]])
			## U Grid
			Grid = lapply(varsTrain,function(x)x$Lt[[i]])
			## Difference by covariate
			for(j in 1:ncovariate){
					Ind_phi = ((j-1)*UnitLen+1) : (j*UnitLen)
					if(L == 1){
							PhiOnGrid = approx(x = LoadingGrid,y = Phi[Ind_phi,1],xout = Grid[[j]] )$y
							Z =  U[[j]] - PhiOnGrid * score[[j]][i,1:L]
						}else{
							PhiOnGrid = apply(Phi[Ind_phi,1:L],2,function(v) approx(x = LoadingGrid,y = v,xout = Grid[[j]] )$y )
							Z =  U[[j]] - PhiOnGrid%*% score[[j]][i,1:L]
						}
					psuLogli[L] = psuLogli[L] + sum(Z^2)/(2*Sigma2[[j]]) + (length(Grid[[j]])/2)*log(pi*2) + (length(Grid[[j]])/2)*log(Sigma2[[j]])
				}
		}
		if(BIC == FALSE){
			psuLogli[L] = psuLogli[L]+ L	
		}else{
			psuLogli[L] = psuLogli[L]+ L*log(nsubj)/2
		}
	}
	Lchoice = which.min(psuLogli)
	Lchoice
}

PredSingle<-function(U,CrMatInfo,EigenInfo){
	ncovariate = length(CrMatInfo$FPCAlist)
	UnitLen = CrMatInfo$FPCAlist[[1]]$workGrid
	for(i in 1:ncovariate){
		PCARes = CrMatInfo$FPCAlist[[i]]
		mu_x = approx(x = PCARes$workGrid, y = PCARes$mu,xout  = U[[i]]$Lt)$y
		Ind = ((i-1)*UnitLen+1) : (i*UnitLen)
		CovMat = CrMatInfo$MultiCrXY[Ind,]
		ginv(CovMat)
	}
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

########################
##function CrWholeMat
##Input: scalar response Y
##       multiple functional process list X
##		 Recorded timePoint tPoint
##		 FPCA settings options
##Output:  Returnlist: a list consists of 
##						FPCAlist: a list of the FPCA results of each X
##						CrMat: The whole Covariance matrix of list X
##						CrMatYZ: a list of the cross covariance matrices of X,Y 
########################
CrWholeMat<-function(Y,X,tPoint,options){
	nsubj = length(X)
	##doing FPCA separately
		for(i in 1:nsubj){
			assign( paste0("PCARes_",i),FPCA(X[[i]],tPoint[[i]],options[[i]]) )
		}
		##choose CrMat size baseon dataType
		TPlength = length( get(paste0("PCARes_",1))$workGrid)
		CrMat = matrix(0,nsubj*TPlength,nsubj*TPlength)
		## fill in the holes
		for (i in 1:nsubj) {
			tmp = ((i-1)*TPlength+1) : (i*TPlength)
			CrMat[tmp,tmp] = get(paste0("PCARes_",i))$smoothedCov
		}
		##calculate the cross cov
		for(i in 1:nsubj){
			for(j in (i):nsubj){
				if(i!=j){
					if(options[[1]]$dataType == "Dense"){
						PCARes_A = get(paste0("PCARes_",i))
						PCARes_B = get(paste0("PCARes_",j))
						CrCovBlock = cov(matrix(unlist(X[[i]]), ncol = length(X[[i]][[1]]), byrow = TRUE),
										matrix(unlist(X[[j]]), ncol = length(X[[j]][[1]]), byrow = TRUE))
						tmp_row = ((i-1)*TPlength+1) : (i*TPlength)
						tmp_col = ((j-1)*TPlength+1) : (j*TPlength)
						CrMat[tmp_row,tmp_col] = CrCovBlock
						CrMat[tmp_col,tmp_row] = t(CrCovBlock)
					}else if(options[[1]]$dataType == "Sparse"){
						PCARes_A = get(paste0("PCARes_",i))
						PCARes_B = get(paste0("PCARes_",j))
						CrCovBlock = GetCrCovYX(bw1 = 0.1, Ly1 = X[[i]], Lt1 = tPoint[[i]],
					 				Ymu1 = approx(x = PCARes_A$workGrid,y = PCARes_A$mu,xout = PCARes_A$workGrid,rule = 2 )$y
					 				, bw2 = 0.1, Ly2 = X[[j]],
				  	 				Lt2 = tPoint[[j]], 
				  	 				Ymu2 =approx(x = PCARes_B$workGrid,y = PCARes_B$mu,xout = PCARes_A$workGrid,rule = 2 )$y)
						tmp_row = ((i-1)*TPlength+1) : (i*TPlength)
						tmp_col = ((j-1)*TPlength+1) : (j*TPlength)
						CrMat[tmp_row,tmp_col] = CrCovBlock$smoothedCC
						CrMat[tmp_col,tmp_row] = t(CrCovBlock$smoothedCC)
					}
				}
			}
		}
		##return FPCA result and covMat
		FPCAlist = list()
		for(i in 1:nsubj){
			FPCAlist = append(FPCAlist,list(get(paste0("PCARes_",i))) )
		}
		CrMatYZ = list()
		for(i in 1:nsubj){
			if(varsOptns[[1]]$dataType == "Dense"){
				##under dense case, directly use matrix calculation
				g_hat = 0
				for(j in 1:n){
		  			g_hat = g_hat+(Y[j]-mean(Y))*(X[[i]][[j]])/n
				}
				CrMatYZ = append(CrMatYZ,list( g_hat ) )
				## warning: if X have missing values involved or not observed on same time point for different curves
				## our method cannot be used here.
			}else if(varsOptns[[1]]$dataType == "Sparse"){
				##under sparse case
				mu_imputed = approx(x = get(paste0("PCARes_",i))$workGrid,y = get(paste0("PCARes_",i))$mu,xout = sort(unique(unlist(tPoint))),rule = 2)$y
				CrCovInfo = GetCrCovYZ(bw = 0.1, Y, Zmu = mean(Y), X[[i]], Lt = tPoint, Ymu = mu_imputed,
							support = get(paste0("PCARes_",i))$workGrid,kern = "gauss")
				CrMatYZ = append(CrMatYZ,list( CrCovInfo$smoothedCC ) )
			}
		}
		Returnlist = list(FPCAlist,CrMat,CrMatYZ)
		names(Returnlist) = c("FPCAlist","MultiCrXY","MultiCrYZ")
			return(Returnlist)	
}

