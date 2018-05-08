#' Function for performing functonal linear regression where the covariates are functions X1(t1),X2(t2),.. and the response is a scalar Y.
#'
#' @param vars A list of input functional covariates with name of "X1", "X2",.. and a scalar response with name "Y". Each field should have two fields: 'Lt', a list (sparse) or a matrix (Dense) specifying the time of observations, and 'Ly', a vetor stand for the observations.
#' @param varsOptns A list of options named by "X1", "X2",... Each filed specify the paramaters that control the corresponding variables. (default: see details of FPCA())
#' @param isNewSub A 1*n vector of 0s or 1s, where n is the total count of subjects. 0 denotes the corresponding subject is only used for estimation and 1 denotes the corresponding subject is only used for prediction. (default: 0's)
#' @param methodSelect The method used for selecting the number of principal components of functional predictors X's used in functional regression , method optns including 'AIC', 'BIC','Basis' and 'FVE'. (default: "FVE") If choose 'FVE' method, parameter 'FVEThreshold' should be implement here. If choose 'Basis' method, 'NumberOfBasis' should be fill in.
#' @param interval A length 2 vector indicate the support grid of the functional covariate, default to be c(0,1)
#'
#' @return A list containing the following fields:
#' \item{BetaList}{A list with fields of estimated beta_XiY(s,t) defiend on [range(Xi),range(Y)]}
#' \item{R2}{Functional R-square.}
#' \item{predictions}{predict value}
#' \item{NOC}{number of component chosen here}
#' \item{EigenV}{FPC eigen value}
#' \item{EigenF}{FPC eigen function (generated from large function)}
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.}
#'
#' @examples 
#' set.seed(1000)
#' #Model: E(Y(t)|X) = int(beta(s,t)*X(s))
#' n <- 200 #number of subjects
#' ngrids <- 51 #number of grids in [0,1] for X(s)
#' ngridt <- 101 #number of grids in [0,1] for Y(t)
#' grids <- seq(0, 1, length.out=ngrids) #regular grids in [0,1] for X(s)
#' gridt <- seq(0, 1, length.out=ngridt) #regular grids in [0,1] for Y(t)
#' #generate X
#' #{1, sqrt(2)*sin(2*pi*s), sqrt(2)*cos(2*pi*t)} are used to generate X.
#' eigenFun <- list(function(s){1 + 0 * s},function(s){sqrt(2) * sin(2*pi*s)},function(s){sqrt(2) * cos(2*pi*s)})
#' sig <- matrix(c(1.5, 0.0, 0.0, 0.9, -.5, 0.1,
#'                 0.0, 1.2, 0.0, -.3, 0.8, 0.4,
#'                 0.0, 0.0, 1.0, 0.4, -.3, 0.7,
#'                 0.9, -.3, 0.4, 2.0, 0.0, 0.0,
#'                 -.5, 0.8, -.3, 0.0, 1.5, 0.0,
#'                 0.1, 0.4, 0.7, 0.0, 0.0, 1.0),
#'                 nrow=6,ncol=6)
#' scoreX <- MASS::mvrnorm(n,mu=rep(0,6),Sigma=sig)
#' scoreX1 <- scoreX[,1:3]
#' scoreX2 <- scoreX[,4:6]
#' basisX1 <- sapply(eigenFun,function(x){x(grids)})
#' latentX1 <- scoreX1 %*% t(basisX1)
#' measErrX1 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
#' denseX1 <- latentX1 + measErrX1
#' basisX2 <- sapply(eigenFun,function(x){x(grids)})
#' latentX2 <- scoreX2 %*% t(basisX2)
#' measErrX2 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
#' denseX2 <- latentX2 + measErrX2
#' #generate Y
#' #beta(s, t) <- sin(2 * pi * s)*cos(2 * pi * t)
#' betaEigen <- list(function(s){1 + 0 * s},function(s){sqrt(2) * sin(2*pi*s)},function(s){sqrt(2) * cos(2*pi*s)})
#' basisY <- c(1,0.1,0.01)
#' latentY <- scoreX1 %*% basisY + scoreX2 %*% basisY
#' measErrY <- sqrt(0.01) * rnorm(n) #0.01 is sigma^2
#' denseY <- latentY + measErrY
#' #======Dense data===============================================
#' timeX <- t(matrix(rep(grids, n),length(grids), n))
#' denseVars <- list(X1 = list(Ly = denseX1, Lt = timeX),
#'                   X2 = list(Ly = denseX2, Lt = timeX),
#'                   Y= denseY)
 
#' resuDense <- FPCRegS(denseVars) 
#' #======Sparse data===============================================
#' sparsity = 5:8
#' sparseX1 <- Sparsify(denseX1, grids, sparsity)
#' sparseX2 <- Sparsify(denseX2, grids, sparsity)
#' sparseVars <- list(X1 = sparseX1, X2 = sparseX2, Y = denseY)
 
#' resuSparse <- FPCRegS(sparseVars, methodSelect=list(method = "FVE",FVEThreshold = 0.9)) #or resuSparse <- FPCReg(vars = sparseVars,varsOptns = list(X1=list(userBwCov = 0.03)))
#' MixedVars <- list(X1 = sparseX1, X2 = list(Ly = denseX2, Lt = timeX), Y = denseY)
#' TestRes3 = FPCRegS(MixedVars)
 
#' basisX2 <- sapply(eigenFun,function(x){x(grids)})
#' latentX2 <- scoreX2 %*% t(basisX2)
#' measErrX2 <- sqrt(0.03) * matrix(rnorm(n * ngrids), n, ngrids) #0.01 is sigma^2.
#' denseX2 <- latentX2 + measErrX2
#'
#' #generate Y
#' #beta(s, t) <- sin(2 * pi * s)*cos(2 * pi * t)
#' betaEigen <- list(function(s){1 + 0 * s},function(s){sqrt(2) * sin(2*pi*s)},function(s){sqrt(2) * cos(2*pi*s)})
#' basisY <- c(1,0.1,0.01)
#' latentY <- scoreX1 %*% basisY + scoreX2 %*% basisY
#' measErrY <- sqrt(0.01) * rnorm(n) #0.01 is sigma^2
#' denseY <- latentY + measErrY
#'
#' #======Dense data===============================================
#' timeX <- t(matrix(rep(grids, n),length(grids), n))
#' denseVars <- list(X1 = list(Ly = denseX1, Lt = timeX),
#'                   X2 = list(Ly = denseX2, Lt = timeX),
#'                   Y= denseY)
#'
#' resuDense <- FPCRegS(denseVars) 
#' #======Sparse data===============================================
#' sparsity = 5:8
#' sparseX1 <- Sparsify(denseX1, grids, sparsity)
#' sparseX2 <- Sparsify(denseX2, grids, sparsity)
#' sparseVars <- list(X1 = sparseX1, X2 = sparseX2, Y = denseY)
#'
#' resuSparse <- FPCRegS(sparseVars, methodSelect=list(method = "FVE",FVEThreshold = 0.9)) #or resuSparse <- FPCReg(vars = sparseVars,varsOptns = list(X1=list(userBwCov = 0.03)))
#' MixedVars <- list(X1 = sparseX1, X2 = list(Ly = denseX2, Lt = timeX), Y = denseY)
#' TestRes3 = FPCRegS(MixedVars)
#' @export


FPCRegS <- function(vars, varsOptns = NULL, isNewSub = NULL, methodSelect = list(method = "FVE",FVEThreshold = 0.95,
	NumberOfBasis = NULL) ,interval = c(0,1)){
	#===============data checking and manipulation===============
	p <- length(vars)-1
	if (p == 0) stop('Too few covariates.')
	if (!'Y' %in% names(vars)|"" %in% names(vars)) stop('Missing name of the response which should be "Y".')
	if (!'X1' %in% names(vars)|"" %in% names(vars)) stop('Missing name of preictors which should be "X1","X2",...')
	if (anyDuplicated(names(vars)) > 0) stop('Duplicated Names.')
	## separate Y and covariates
	if ('Y' %in% names(vars)) {
		Y = vars['Y']
		vars <- vars[names(vars) != 'Y'] 
	}
	#Set Lt in [0,1] for dense data without Lt. 
	for (i in 1:(p)) {
		if (is.null(vars[[i]]$Lt) & is.matrix(vars[[i]]$Ly)) {
			warning('time variate missing, choose the equally spaced grid here.')
			n1 <- dim(vars[[i]]$Ly)[1]
			n2 <- dim(vars[[i]]$Ly)[2]
			vars[[i]]$Lt <- t(matrix(rep(seq(0, 1, length.out = n2), n1), n2, n1))
		}
	}
	#set options 
	#varsOptnsDef is the default
	#specific default setting for optns now and output grids for sparse is fixed at 51 now due to the cross-cov.
	if (is.list(varsOptns) | is.null(varsOptns)) {
		for(i in 1:p){
				Dense = ifelse(is.matrix(vars[[i]]$Ly)&is.matrix(vars[[i]]$Lt),1,0)
				## For the case with regular design, also think as dense
				if( sum(unlist(lapply(vars[[i]]$Lt,function(x){ if(length(x) == length(vars[[i]]$Lt[[1]])) sum(x - vars[[i]]$Lt[[1]]) }))) == 0  ){Dense = 1}
				if (Dense == 1) {optns <- list(dataType = "Dense",error = TRUE, kernel='gauss', nRegGrid=51, useBinnedData='OFF')} else {optns <- list(dataType = "Sparse", error = TRUE, kernel = 'gauss' ,nRegGrid = 51, useBinnedData = 'OFF',methodBwCov = 'GMeanAndGCV')} 	
				if (i == 1) {varsOptnsDef <- list(optns)} else {varsOptnsDef <- c(varsOptnsDef, list(optns))}
		}
		names(varsOptnsDef) <- names(vars)
		if(is.null(varsOptns)){varsOptns <- varsOptnsDef
		}} else if (!is.list(varsOptns)) {stop('wrong input of varsOptns.')}
	if (!sum(names(varsOptns) %in% names(vars)) == length(names(varsOptns))) stop('Check names of varsOptns which should be X1,X2..Y.')
	for (i in 1:(p)) {
		Dense = ifelse(is.matrix(vars[[i]]$Ly)&is.matrix(vars[[i]]$Lt),1,0)
		name1=names(varsOptns)
		if (!names(vars)[i]%in%names(varsOptns)) {varsOptns <- c(varsOptns, list(optns)); names(varsOptns) = c(name1, names(vars)[i])}
		if (is.null(varsOptns[names(vars)[i]]$dataType) & Dense == 0) {varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], dataType = "Sparse")}else if(is.null(varsOptns[names(vars)[i]]$dataType) & Dense==1){varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], dataType = "Dense")}
		if (is.null(varsOptns[names(vars)[i]]$nRegGrid) & Dense == 0) {varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], nRegGrid = 51)}
		if (is.null(varsOptns[names(vars)[i]]$error)) {varsOptns[[names(vars)[i]]] <- c(varsOptns[[names(vars)[i]]], error = TRUE)}
		if (Dense == 1) {varsOptns[[names(vars)[i]]]$nRegGrid <- ncol(vars[[i]]$Ly)}
		}
	varsOptns <- varsOptns[names(vars)]

	#nRegGrids determine the numbers of ouput grids 
	#if (Dense == 0) {nRegGrids <- sapply(varsOptns, function(x){x$nRegGrid})}
	for ( i in 1:(p)) {
		if (!'Ly' %in% names(vars[[i]])) stop('Insert the name "Ly" for the predictors and response to indicate time of data.')
		if (!'Lt' %in% names(vars[[i]])) stop('Insert the name "Lt" for the predictors and response to indicate time of data.')
      	if (! length(names(vars[[i]])) == 2) stop('Check data.')
		}
	for (i in 1:(p)) {
		if (is.list(vars[[i]]$Lt) == 0) {vars[[i]]$Lt <- lapply(1:nrow(vars[[i]]$Lt), function(j) vars[[i]]$Lt[j, ])}
		if (is.list(vars[[i]]$Ly) == 0) {vars[[i]]$Ly <- lapply(1:nrow(vars[[i]]$Ly), function(j) vars[[i]]$Ly[j, ])}
	}
	#list dense data for using HandleNumericsAndNAN and demeanFuc func.
	for (i in 1:(p)) {
		if (is.list(vars[[i]]$Lt) == 0) {vars[[i]]$Lt <- lapply(1:nrow(vars[[i]]$Lt), function(j) vars[[i]]$Lt[j, ])}
		if (is.list(vars[[i]]$Ly) == 0) {vars[[i]]$Ly <- lapply(1:nrow(vars[[i]]$Ly), function(j) vars[[i]]$Ly[j, ])}
		}
	vars[sapply(vars, is.list)]<-lapply(vars[sapply(vars, is.list)], function(v) HandleNumericsAndNAN(v[['Ly']], v[['Lt']]))
	if (is.null(isNewSub)) {isNewSub <- rep(0, length(vars[[1]]$Lt))}
	varsTrain <- vars
	for (i in 1:(p)) {
		varsTrain[[i]]$Lt <- vars[[i]]$Lt[which(isNewSub == 0)]
		varsTrain[[i]]$Ly <- vars[[i]]$Ly[which(isNewSub == 0)]
		}
	Y_Train = unlist(Y)[which(isNewSub == 0)]
	#===============population parameters===============	
	demeanedRes <- demeanFuc(p-1, varsTrain, kern='gauss', varsOptns) #Centered predictors. Using gauss for demeanFuc, but may be relaxed.
	varsTrain <- demeanedRes[['xList']]
	muList <- demeanedRes[['muList']]
	#### make list
	TPlist<-NULL
	XList<-NULL
	for(i in 1:p){
		if(i==1){
			TPlist = list(varsTrain$X1$Lt)
			XList = list(varsTrain$X1$Ly)
			}else{
				TPlist = append(TPlist,list(varsTrain[[i]]$Lt))
				XList = append(XList,list(varsTrain[[i]]$Ly))
			}
	}
	GetMatrixInfo = CrWholeMat(Y_Train,XList,TPlist,varsOptns)
	TimePoint = GetMatrixInfo$FPCAlist[[1]]$workGrid
	Eigen_Decompose = eigen(GetMatrixInfo$MultiCrXY * (TimePoint[2]-TimePoint[1]))
	Eigen_Value = Eigen_Decompose$values[Eigen_Decompose$values>0]
	Eigen_fun = Eigen_Decompose$vectors[,Eigen_Decompose$values>0]
	AdjForInt = apply(Eigen_fun,2,function(x){MultiNorm(x^2,p,TimePoint,interval)})
	#Eigen_Value = Eigen_Value*AdjForInt
	 for(i in 1:length(AdjForInt)){
	 	Eigen_fun[,i] = Eigen_fun[,i]/sqrt(TimePoint[2]-TimePoint[1])
	 }
		if(methodSelect$method == "FVE"){
			if(is.null(methodSelect$FVEThreshold)){stop("FVEThreshold is needed for method FVE")}
			Ratio = cumsum(Eigen_Value)/sum(Eigen_Value)
			L = sum(Ratio < methodSelect$FVEThreshold)+1
		}else if(methodSelect$method == "Basis"){
			L = methodSelect$NumberOfBasis
		}else if(methodSelect$method == "AIC"){
			L = AICchoice(varsTrain,GetMatrixInfo,Eigen_Decompose)
		}else if(methodSelect$method == "BIC"){
			L = AICchoice(varsTrain,GetMatrixInfo,Eigen_Decompose,BIC = TRUE)
		}else{
			stop('no such method select method implemented.')
		}
		beta = 0
		score = rep(0,L)
		for(i in 1:L){
			nominater = MultiNorm(unlist(GetMatrixInfo$MultiCrYZ) * Eigen_fun[,i] , p, TimePoint,interval)
			score[i] = nominater/Eigen_Value[i]
			beta = beta + score[i]*Eigen_fun[,i]
		}
		###Split beta
		Beta = list()
		for(i in 1:p){
			l = length(Eigen_Decompose$vectors[,i])/p
			Ind_beta = ((i-1)*l+1) : (i*l)
			Beta = append(Beta,list(beta[Ind_beta]))
		}
	Betalist = list(list(Beta),list(TimePoint),list(score))
	names(Betalist) = c("estiBeta","TimePoint","b_score")
	#####R2 score
	R2 = sum( score^2*Eigen_Value[1:L] ) / (sum( score^2*Eigen_Value[1:L] )+ sum(unlist(lapply(GetMatrixInfo$FPCAlist,function(x){x$sigma2}))) )

	#####predictions
	###so far we use int(X1beta1)+int(X2beta2)+... 
	if(sum(isNewSub) != 0){
		varsPred <- vars
		Y_Pred = mean(Y_Train)
		Dense = 1
		for(i in 1:p){
			varsPred[[i]]$Lt <- vars[[i]]$Lt[which(isNewSub == 1)]
			varsPred[[i]]$Ly <- vars[[i]]$Ly[which(isNewSub == 1)]
			DenseI = ifelse(is.matrix(varsPred[[i]]$Ly)&is.matrix(varsPred[[i]]$Lt),1,0)
			## For the case with regular design, also think as dense
			if( sum(unlist(lapply(varsPred[[i]]$Lt,function(x){ if(length(x) == length(varsPred[[i]]$Lt[[1]])) sum(x - varsPred[[i]]$Lt[[1]]) }))) == 0  ){DenseI = 1}
			Dense = Dense & DenseI
		}
		if( Dense == 1){
			for(i in 1:p){
				Y_Pred = Y_Pred + unlist(lapply(varsPred[[i]]$Ly, y<-function(y) {FakeInt(Beta[[i]]*approx(x = varsPred[[i]]$Lt[[1]],y = y,xout = GetMatrixInfo$FPCAlist[[i]]$workGrid ,rule = 2)$y,TimePoint,interval)}))	
			}}else{
					Kt = min(floor(L/2),length(GetMatrixInfo$FPCAlist[[i]]$workGrid))
					#PredScore = predict.FPCA(GetMatrixInfo$FPCAlist[[i]],newLy = varsPred[[i]]$Ly,newLt = varsPred[[i]]$Lt,xiMethod = 'CE',K = Kt)
					#X_Imp = GetMatrixInfo$FPCAlist[[i]]$phi[,1:Kt]%*%t(PredScore)
					CEscore = matrix(0,length(varsPred[[1]]$Ly),length(Eigen_Value))
					for(j in 1:length(varsPred[[1]]$Ly)){
					  CEscore[j,] = GetCE_Mul(v = lapply(varsPred,function(x){x$Ly[[j]]}) ,tp = lapply(varsPred,function(x){x$Lt[[j]]}),MatrixInfo = GetMatrixInfo,Lambda = Eigen_Value,Phi = Eigen_fun)
					}
					#Y_Pred = Y_Pred + apply( X_Imp ,2,function(x) FakeInt(x*Beta[[i]], GetMatrixInfo$FPCAlist[[i]]$workGrid,interval)   )
					Y_Pred = Y_Pred + apply(CEscore,1,function(x) sum(x[1:L]*score[1:L]) ) 
				}
		returnList = list(Betalist,R2,Y_Pred,L,Eigen_Value[1:L],Eigen_fun[,1:L])
	names(returnList) = c("Betalist","R2","predictions","NOC","EigenV","EigenF")
	returnList		
	}else{
		Y_Pred = NULL
	}
		returnList = list(Betalist,R2,Y_Pred,L,Eigen_Value)
		names(returnList) = c("Betalist","R2","predictions","NOC","Eigen")
		returnList
}


