#library(Matrix)
#library(fields)#interp.surface.grid

FPCReg <- function(vars,m_optns=NULL,isNewSub=NULL,method='AIC',FVEthreshold=0.98,alpha=0.05){
	#===============data checking and manipulation
	orinames=names(vars)
	p<-length(vars)-1
	if(p==0) stop('Too few covariates.')
	if(!'Y' %in% names(vars)|"" %in% names(vars)) stop('Missing name of the response which should be "Y".')
	if(!'X1' %in% names(vars)|"" %in% names(vars)) stop('Missing name of preictors which should be "X1","X2",...')
	if(anyDuplicated(names(vars))>0) stop('Duplicated Names.')
	if('Y' %in% names(vars)) {vars<-c(vars[names(vars) != 'Y'],vars['Y'])}else if(names(vars)[length(vars)]== ''){names(vars)[length(vars)]<-'Y'} #seems unnecessary
	#names(vars)<-c(paste0('X',seq_len(length(vars)-1)),'Y')
	
	for(i in 1:(p+1)){
		if(is.null(vars[[i]]$Lt)&is.matrix(vars[[i]]$Ly)){
			n1=dim(vars[[i]]$Ly)[1]
			n2=dim(vars[[i]]$Ly)[2]
			vars[[i]]$Lt=t(matrix(rep(seq(0,1,length.out=n2),n1),n2,n1))
			}
		}
	
	if(sum(sapply(vars,function(x){is.matrix(x$Ly)})*sapply(vars,function(x){is.matrix(x$Lt)}))){Dense=1}else{Dense=0}
	if(Dense==1){optns=list(dataType="Dense",error=1,kernel='gauss',nRegGrid=51,useBinnedData='OFF')}else{optns=list(dataType="Sparse",error=TRUE,kernel='gauss',nRegGrid=51,useBinnedData='OFF')}
	#optns=list(dataType="Sparse",error=TRUE,kernel='gauss',nRegGrid=51,useBinnedData='OFF')
	for(i in 1:(p+1)){
		if(i==1){m_optns_temp=list(optns)}else{m_optns_temp=c(m_optns_temp,list(optns))}
		}
	names(m_optns_temp)=names(vars)
	#if(m_optns=='default'){m_optns=m_optns_temp}else if(!is.list(m_optns)){stop('wrong inpiut of m_optns.')}
	if(is.null(m_optns)){m_optns=m_optns_temp}else if(!is.list(m_optns)){stop('wrong inpiut of m_optns.')}
	if(!sum(names(m_optns) %in% names(vars))==length(names(m_optns))) stop('Check names of m_optns which should be X1,X2..Y.')
	for(i in 1:(p+1)){
		name1=names(m_optns)
		if(!names(vars)[i]%in%names(m_optns)){m_optns=c(m_optns,list(optns));names(m_optns)=c(name1,names(vars)[i])}
		if(is.null(m_optns[names(vars)[i]]$dataType)&Dense==0){m_optns[[names(vars)[i]]]=c(m_optns[[names(vars)[i]]],dataType = "Sparse")}else if(is.null(m_optns[names(vars)[i]]$dataType)&Dense==1){m_optns[[names(vars)[i]]]=c(m_optns[[names(vars)[i]]],dataType = "Dense")}
		if(is.null(m_optns[names(vars)[i]]$nRegGrid)&Dense==0){m_optns[[names(vars)[i]]]=c(m_optns[[names(vars)[i]]],nRegGrid = 51)}
		if(Dense==1){m_optns[[names(vars)[i]]]$nRegGrid=dim(vars[[i]]$Ly)[2]}
		}
	m_optns=m_optns[names(vars)]
	
	#if(sum(sapply(m_optns,function(x){!x$dataType=="Dense"}))==0){Dense=1}else{Dense=0}
	if(Dense==0){nRegGrids=sapply(m_optns,function(x){x$nRegGrid})}

	for( i in 1:(p+1)){
		if(!'Ly' %in% names(vars[[i]])) stop('Insert the name "Ly" for the predictors and response to indicate time of data.')
		if(!'Lt' %in% names(vars[[i]])) stop('Insert the name "Lt" for the predictors and response to indicate time of data.')
      	if(! length(names(vars[[i]]))==2) stop('Check data.')
		}

	for(i in 1:(p+1)){
		if(is.list(vars[[i]]$Lt)==0){vars[[i]]$Lt=lapply(1:nrow(vars[[i]]$Lt), function(j) vars[[i]]$Lt[j,])}
		if(is.list(vars[[i]]$Ly)==0){vars[[i]]$Ly=lapply(1:nrow(vars[[i]]$Ly), function(j) vars[[i]]$Ly[j,])}
		} #list data for HandleNumericsAndNAN and demean func.

	vars[sapply(vars,is.list)]<-lapply(vars[sapply(vars, is.list)], function(v) HandleNumericsAndNAN(v[['Ly']], v[['Lt']]))

	if(is.null(isNewSub)){isNewSub=rep(1,length(vars[[1]]$Lt))}
	vars_train=vars
	vars_test=vars
	for(i in 1:(p+1)){
		vars_train[[i]]$Lt=vars[[i]]$Lt[which(isNewSub==1)]
		vars_train[[i]]$Ly=vars[[i]]$Ly[which(isNewSub==1)]	
		#vars_test[[i]]$Lt=vars[[i]]$Lt[which(isNewSub==0)]
		#vars_test[[i]]$Ly=vars[[i]]$Ly[which(isNewSub==0)]
		}

	#===============population parameters
	intLen=array(0,p+1);grid_n=array(0,p+1)
	for(i in 1:(p+1)){
		intLen[i]=max(c(unlist(vars_train[[i]]$Lt)))-min(c(unlist(vars_train[[i]]$Lt)))	
		grid_n[i]=if(Dense==1){length(vars_train[[i]]$Ly[[1]])}else{nRegGrids[i]}
		}
	intLen_x=intLen[1:p]		  #time window lengthes
	grid_x=grid_n[1:p]            #number of grids for predictors
	brk1=c(0,cumsum(grid_x))      #break points for matrix products
	mx_wet=wet_matrix(p,intLen_x,grid_x,brk1) #library(Matrix) for bdiag()

	demeanedRes=demean(p,vars_train,kern='gauss',m_optns) #just use gauss for demean, ma be relax
	vars_train<-demeanedRes[['xList']]
	muList<-demeanedRes[['muList']]
	if(Dense==1){
		for(i in 1:(p+1)){
			vars_train[[i]]$Lt=do.call(rbind, vars_train[[i]]$Lt)
			vars_train[[i]]$Ly=do.call(rbind, vars_train[[i]]$Ly)
			}
		}

	mvCov=if(Dense==1){dense_Cov(p,vars_train,brk1,mx_wet,grid_n,m_optns)}else{sparse_Cov(p,vars_train,brk1,m_optns,nRegGrids,vars,muList,grid_x,mx_wet)}
	preCov=mvCov[['preCov']]
	croCov=mvCov[['croCov']]
	yCov=mvCov[['yCov']]
	sigma=mvCov[['sigma']]
	n_positive=mvCov[['n_positive']]
	evalue=mvCov[['evalue']]	
	evector=mvCov[['evector']]
	
	if(Dense==0){
		varsDemean=mvCov[['varsDemean']]
		subMx=mvCov[['subMx']]	
		meanVars=mvCov[['meanVars']]	
		phimx=mvCov[['phimx']]	
		meanX=mvCov[['meanX']]	
		diagmx=mvCov[['diagmx']]	
		n_positive=mvCov[['n_positive']]
		inv_mx=lapply(subMx,function(x){as.matrix(solve(x))})
		productrel_1=mapply(function(X,Y){X%*%Y},X=inv_mx,Y=varsDemean)
		productrel_2=mapply(function(X,Y){diag(evalue[1:n_positive])%*%X%*%Y},X=phimx,Y=productrel_1,SIMPLIFY=FALSE)
	}else{
		sigmaVar=mvCov[['sigmaVar']]
		meanVars=list()
		for(i in 1:p){
			gr_x=seq(min(unlist(vars[[i]]$Lt)),max(unlist(vars[[i]]$Lt)),length.out=grid_x[i])
			if(i==1){MxAll_0=t(sapply(vars[[1]]$Ly,c))}else{MxAll_0=cbind(MxAll_0,t(sapply(vars[[i]]$Ly,c)))}
			if(i==1){muListX=muList[[1]](gr_x)}else{muListX=c(muListX,muList[[i]](gr_x))}
			meanVars[[i]]=muList[[i]](gr_x)
			}
		MxAll=apply(MxAll_0,1,'-',muListX)
		MxAll=t(MxAll)
		Mxscore=mat_prod(MxAll,evector,brk1,intLen_x)
		}

	if(Dense==0){n_pcomp=sparse_pcacomp(productrel_2,varsDemean,method,FVEthreshold,evalue,n_positive,isNewSub,diagmx,phimx)}else{n_pcomp=dense_pcacomp(method,MxAll,Mxscore,evector,evalue,isNewSub,sigmaVar,n_positive,FVEthreshold)}
	es_beta=array(0,c(dim(croCov)))
	for(i in 1:n_pcomp){
		es_beta=es_beta+(1/evalue[i])*evector[,i]%*%mat_prod(t(evector[,i]),croCov,brk1,intLen)
		}
	esti_beta=list()
	for(i in 1:p){
		esti_beta[[i]]=es_beta[(brk1[i]+1):brk1[i+1],]
		}
	names(esti_beta)=lapply(names(vars)[1:p],function(x){paste("beta_",x,"_Y",sep='')})

	intbetaphi=mat_prod(t(es_beta),as.matrix(evector[,1:n_pcomp]),brk1,intLen)
	diagYcov=diag(yCov)
	varY=mat_prod(t(diagYcov),as.matrix(rep(1,grid_n[p+1])),c(0,grid_n[p+1]),intLen[p+1])
	varYX=sum(mat_prod(t(intbetaphi^2),as.matrix(rep(1,grid_n[p+1])),c(0,grid_n[p+1]),intLen[p+1])*evalue[1:n_pcomp])
	Q=min(varYX/varY,1)

	r2=apply(intbetaphi^2%*%diag(evalue[1:n_pcomp]),1,sum)/diagYcov
	r2=sapply(r2,function(x){min(x,1)})
	R2=mat_prod(t(r2),as.matrix(rep(1,grid_n[p+1])),c(0,grid_n[p+1]),intLen[p+1])/intLen[p+1]

	#===============prediction for Dense and Sparse
	if(Dense==1){
		gr_y=seq(min(unlist(vars[[p+1]]$Lt)),max(unlist(vars[[p+1]]$Lt)),length.out=grid_x[1])
		meanY=muList[[p+1]](gr_y)
		meanVars[[p+1]]=muList[[p+1]](gr_y);names(meanVars)=names(vars)
		predictY=meanY+intbetaphi%*%t(Mxscore[,1:n_pcomp])
		predictY=t(predictY)
		emx=diag(evalue[1:n_pcomp])
		for(i in 1:p){
			if(i==1){
				sigmaVar=c(rep(sigma[1],grid_n[1]))}else{sigmaVar=c(sigmaVar,rep(sigma[i],grid_n[i]))}
			}
		SigmaY=preCov+diag(sigmaVar)
		Omega=emx-emx%*%t(evector[,1:n_pcomp])%*%solve(SigmaY)%*%evector[,1:n_pcomp]%*%t(emx)
		pwvar=diag(intbetaphi%*%Omega%*%t(intbetaphi))
		upCI=t(apply(predictY,1,'+',qnorm(1-alpha/2)*sqrt(pwvar)))
		lwCI=t(apply(predictY,1,'-',qnorm(1-alpha/2)*sqrt(pwvar)))
		CI_Y=list(upCI=upCI,lwCI=lwCI)
	}else{
		score=lapply(productrel_2,function(x){x[1:n_pcomp]})		
		imputeY=lapply(score,function(x){intbetaphi%*%x})
		gr_y=seq(min(unlist(vars[[p+1]]$Lt)),max(unlist(vars[[p+1]]$Lt)),length.out=grid_x[1])
		meanY=muList[[p+1]](gr_y)
		meanVars[[p+1]]=muList[[p+1]](gr_y);names(meanVars)=names(vars)
		predmean=mat_prod(t(es_beta),as.matrix(meanX),brk1,intLen_x)
		predictY=lapply(imputeY,function(x){x+meanY})
		emx=diag(evalue[1:n_pcomp])
		H=lapply(phimx,function(x){diag(evalue[1:n_pcomp])%*%x[1:n_pcomp,]})
		Omega=mapply(function(X,Y){emx-X%*%Y%*%t(X)},X=H,Y=inv_mx,SIMPLIFY=FALSE)
		pwvar=lapply(Omega,function(x){diag(intbetaphi%*%x%*%t(intbetaphi))})
		upCI=mapply(function(X,Y){X+qnorm(1-alpha/2)*sqrt(Y)},X=predictY,Y=pwvar)
		lwCI=mapply(function(X,Y){X-qnorm(1-alpha/2)*sqrt(Y)},X=predictY,Y=pwvar)
		CI_Y=lapply(1:dim(lwCI)[2],function(x){CI=cbind(lwCI[,x],upCI[,x]);colnames(CI)=c('lwCI','upCI');return(CI)})
		}

	res<-list(esti_beta=esti_beta,predictY=predictY,CI_Y=CI_Y,Q=Q,R2=R2,sigma=sigma,meanVars=meanVars,n_pcomp=n_pcomp)
	res
	}

demean <- function(p,vars_train,kern,m_optns) {
	for(i in 1:(p+1)){
	      userBwMu_temp=SetOptions(vars_train[[i]]$Ly,vars_train[[i]]$Lt,m_optns[[i]])$userBwMu
		vars_train[[i]]=list(Lt=vars_train[[i]]$Lt,Ly=vars_train[[i]]$Ly,userBwMu=userBwMu_temp)
		}
	tmp <- lapply(vars_train, function(x) {
		Tin <- sort(unique(unlist(x[['Lt']])))
		xmu <- GetSmoothedMeanCurve(x[['Ly']],x[['Lt']],Tin,Tin[1],list(userBwMu=x[['userBwMu']], kernel=kern))[['mu']]
		muFun<-approxfun(Tin,xmu)
		x[['Ly']] <- lapply(1:length(x[['Ly']]),function(i) x[['Ly']][[i]]- muFun(x[['Lt']][[i]]))
		xmu<-muFun
    		list(x=x, mu=xmu)
  		})
  	xList <- lapply(tmp, `[[`, 'x')
  	muList <- lapply(tmp, `[[`, 'mu')
  	list(xList = xList, muList = muList)
	}

mat_prod<-function(A,B,gridbrk,intlen){      
	brk1=gridbrk
	oupmat=array(0,c(dim(A)[1],dim(B)[2]))
	for(i in 1:(length(brk1)-1)){
		a1=brk1[i]+1
		a2=brk1[i+1]
		A[,a1]=A[,a1]/2;A[,a2]=A[,a2]/2
		oupmat=oupmat+A[,a1:a2]%*%B[a1:a2,]*intlen[i]/(a2-a1)
		}
	return(oupmat)
	}

wet_matrix=function(p,intLen_x,grid_x,brk1){    
	for(i in 1:p){
		if(i==1){mx_wet=diag(intLen_x[1]/(grid_x[1]-1),grid_x[1])}else{mx_wet=bdiag(mx_wet,diag(intLen_x[i]/(grid_x[i]-1),grid_x[i]))}
		mx_wet[brk1[i]+1,brk1[i]+1]=mx_wet[brk1[i]+1,brk1[i]+1]/2
		mx_wet[brk1[i+1],brk1[i+1]]=mx_wet[brk1[i+1],brk1[i+1]]/2
		}
	mx_wet=sqrt(as.matrix(mx_wet))
	return(mx_wet)
	}

dense_Cov=function(p,vars_train,brk1,mx_wet,grid_n,m_optns){
	indRaw=vars_train[[1]]$Ly
	if(m_optns[[1]]$error==1){covDense=GetCovDense(indRaw,optns=m_optns[[1]]);mx_temp=covDense$smoothCov;sigma=covDense$sigma2}else{mx_temp=cov(indRaw)}
	if(p>=2){
		for(i in 2:p){	
			indRaw=vars_train[[i]]$Ly		
			if(m_optns[[i]]$error==1){covDense=GetCovDense(indRaw,optns=m_optns[[i]]);mx_temp=bdiag(mx_temp,covDense$smoothCov);sigma=c(sigma,covDense$sigma2)}else{mx_temp=bdiag(mx_temp,cov(indRaw))}	
			}

		for(i in 1:p){
			grid_brk1=brk1[i]+1
			grid_brk2=brk1[i+1]
			for(j in 1:p){
				if(j>i){
					grid_brk3=brk1[j]+1
					grid_brk4=brk1[j+1]
					mx_temp[grid_brk1:grid_brk2,grid_brk3:grid_brk4]=cov(vars_train[[i]]$Ly,vars_train[[j]]$Ly)
					mx_temp[grid_brk3:grid_brk4,grid_brk1:grid_brk2]=t(as.matrix(mx_temp[grid_brk1:grid_brk2,grid_brk3:grid_brk4]))
					}
				}
			}
		}
	preCov=as.matrix(mx_temp)
	spec_preCov=eigen(mx_wet%*%preCov%*%mx_wet)
	evalue=spec_preCov$value
	evalue=sapply(evalue,function(x){max(x,0)})
	n_positive=sum(evalue>0)
	evector=diag(diag(1/mx_wet))%*%spec_preCov$vector
	preCov=evector%*%diag(evalue)%*%t(evector)

	for(i in 1:p){
		if(i==1){croCov=cov(vars_train[[i]]$Ly,vars_train[[p+1]]$Ly)}else{croCov=rbind(croCov,cov(vars_train[[i]]$Ly,vars_train[[p+1]]$Ly))}
		}
	indRaw=vars_train[[p+1]]$Ly
	if(m_optns[[p+1]]$error==1){yCov=GetCovDense(indRaw,optns=m_optns[[p+1]])$smoothCov}else{yCov=cov(indRaw)}
	for(i in 1:p){
		if(i==1){sigmaVar=c(rep(sigma[1],grid_n[1]))}else{sigmaVar=c(sigmaVar,rep(sigma[i],grid_n[i]))}
		}
	return(list(preCov=preCov,croCov=croCov,yCov=yCov,sigma=sigma,evalue=evalue,evector=evector,sigmaVar=sigmaVar,n_positive=n_positive))
	}

sparse_Cov=function(p,vars_train,brk1,m_optns,nRegGrids,vars,muList,grid_x,mx_wet){
	Ly=vars_train[[1]]$Ly;Lt=vars_train[[1]]$Lt;
	optns=SetOptions(Ly,Lt,m_optns[[1]])
	obsGrid=sort(unique(c(unlist(Lt))))
	meanCr=array(0,c(length(obsGrid)))
	meanCrAll=list(meanCr)
	rfpca=FPCA(Ly,Lt,optns)
	mx_temp=rfpca$smoothedCov
	sigma=rfpca$sigma2
	m_bwMu=list(rfpca$bwMu)
	m_bwCov=list(rfpca$bwCov)
	#persp(1:51,1:51,mx_temp1,phi = 45, theta = 45)
	#persp(1:51,1:51,preCov[52:102,52:102],phi = 45, theta = 45)

	if(p>=2){
		for(i in 2:p){
			Ly=vars_train[[i]]$Ly;Lt=vars_train[[i]]$Lt;
			optns=SetOptions(Ly,Lt,m_optns[[i]])
			obsGrid=sort(unique(c(unlist(Lt))))
			meanCr=array(0,c(length(obsGrid)))
			meanCrAll[[i]]=meanCr
			rfpca=FPCA(Ly,Lt,optns)
			mx_temp1=rfpca$smoothedCov
			assign(paste("sigma_",i,sep=""),rfpca$sigma2)
			sigma=c(sigma,get(paste("sigma_",i,sep="")))
			mx_temp=bdiag(mx_temp,mx_temp1)	
			m_bwMu[[i]]=rfpca$bwMu
			m_bwCov[[i]]=rfpca$bwCov
			}

		for(i in 1:p){
			grid_brk1=brk1[i]+1
			grid_brk2=brk1[i+1]
			for(j in 1:p){
				if(j>i){
					grid_brk3=brk1[j]+1
					grid_brk4=brk1[j+1]
					mx_temp[grid_brk1:grid_brk2,grid_brk3:grid_brk4]=GetCrCovYX(bw1=m_bwCov[[i]],bw2=m_bwCov[[j]],Ly1=vars_train[[i]]$Ly,Lt1=vars_train[[i]]$Lt,Ymu1=meanCrAll[[i]],Ly2=vars_train[[j]]$Ly,Lt2=vars_train[[j]]$Lt,Ymu2=meanCrAll[[j]])$smoothedCC
					mx_temp[grid_brk3:grid_brk4,grid_brk1:grid_brk2]=t(as.matrix(mx_temp[grid_brk1:grid_brk2,grid_brk3:grid_brk4]))
					}
				}
			}
		}
	preCov=as.matrix(mx_temp)
	spec_preCov=eigen(mx_wet%*%preCov%*%mx_wet)
	evalue=spec_preCov$value
	evalue=sapply(evalue,function(x){max(x,0)})
	n_positive=sum(evalue>0)
	evector=diag(diag(1/mx_wet))%*%spec_preCov$vector
	preCov=evector%*%diag(evalue)%*%t(evector)
	
	Ly=vars_train[[p+1]]$Ly;Lt=vars_train[[p+1]]$Lt;
	optns=SetOptions(Ly,Lt,m_optns[[p+1]])
	obsGrid=sort(unique(c(unlist(Lt))))
	yfpca=FPCA(Ly,Lt,optns)
	yCov=yfpca$fittedCov
	m_bwMu[[p+1]]=yfpca$bwMu
	m_bwCov[[p+1]]=yfpca$bwCov
	meanCrY=array(0,c(length(obsGrid)))
	for(i in 1:p){
		if(i==1){croCov=GetCrCovYX(bw1=m_bwCov[[i]],bw2=m_bwCov[[p+1]],Ly1=vars_train[[i]]$Ly,Lt1=vars_train[[i]]$Lt,Ymu1=meanCrAll[[i]],Ly2=vars_train[[p+1]]$Ly,Lt2=vars_train[[p+1]]$Lt,Ymu2=meanCrY)$smoothedCC}
		else{croCov=rbind(croCov,GetCrCovYX(bw1=m_bwCov[[i]],bw2=m_bwCov[[p+1]],Ly1=vars_train[[i]]$Ly,Lt1=vars_train[[i]]$Lt,Ymu1=meanCrAll[[i]],Ly2=vars_train[[p+1]]$Ly,Lt2=vars_train[[p+1]]$Lt,Ymu2=meanCrY)$smoothedCC)}
		}
	#==========================
	for(i in 1:p){
		varsDemean0=mapply('-',vars[[i]]$Ly,lapply(vars[[i]]$Lt,muList[[i]]),SIMPLIFY=FALSE)
		if(i==1){varsDemean=varsDemean0}else{varsDemean=mapply(function(X,Y){c(X,Y)},X=varsDemean,Y=varsDemean0)}
		}
	pnumber0=rep(list(0),lengthVars(vars))
	for(i in 1:p){
		pnumberr=sapply(vars[[i]]$Lt,length)
		pnumber0=mapply(c,pnumber0,pnumberr, SIMPLIFY=FALSE)
		}
	pnumber=lapply(pnumber0,cumsum)
	if(p==1){diagmx=lapply(pnumber0,function(x){C=diag(rep(sigma[1],x[2]))})}else{diagmx=lapply(pnumber0,function(x){C=diag(rep(sigma[1],x[2]));for(i in 2:p){C=bdiag(C,diag(rep(sigma[i],x[i+1])))};return(as.matrix(C))})}

	a1=brk1[1]+1
	a2=brk1[1+1]
	gr_i=seq(min(unlist(vars[[1]]$Lt)),max(unlist(vars[[1]]$Lt)),length.out=grid_x[1])
	obj<-list(x=gr_i,y=gr_i,z=preCov[a1:a2,a1:a2])
	subMx=mapply(function(X,Y){interp.surface.grid(obj,list(x=X,y=Y))},X=vars[[1]]$Lt,Y=vars[[1]]$Lt)[3,]
	assign(paste("subMx_",i,i,sep=""),subMx)
	if(p>=2){
		for(i in 2:p){
			a1=brk1[i]+1
			a2=brk1[i+1]
			gr_i=seq(min(unlist(vars[[i]]$Lt)),max(unlist(vars[[i]]$Lt)),length.out=grid_x[i])
			obj<-list(x=gr_i,y=gr_i,z=preCov[a1:a2,a1:a2])
			assign(paste("subMx_",i,i,sep=""),mapply(function(X,Y){interp.surface.grid(obj,list(x=X,y=Y))},X=vars[[i]]$Lt,Y=vars[[i]]$Lt)[3,])
			subMx=mapply(function(X,Y){bdiag(X,Y)},subMx,get(paste("subMx_",i,i,sep="")))	
			}
	
		for(i in 1:p){
			grid_brk1=brk1[i]+1
			grid_brk2=brk1[i+1]
			for(j in 1:p){
				if(j>i){
					grid_brk3=brk1[j]+1
					grid_brk4=brk1[j+1]
					gr_i=seq(min(unlist(vars[[i]]$Lt)),max(unlist(vars[[i]]$Lt)),length.out=grid_x[i])
					gr_j=seq(min(unlist(vars[[j]]$Lt)),max(unlist(vars[[j]]$Lt)),length.out=grid_x[j])
					obj<-list(x=gr_i,y=gr_j,z=preCov[grid_brk1:grid_brk2,grid_brk3:grid_brk4])
					assign(paste("subMx_",i,j,sep=""),mapply(function(X,Y){interp.surface.grid(obj,list(x=X,y=Y))},X=vars[[i]]$Lt,Y=vars[[j]]$Lt)[3,])
					subMxx=get(paste("subMx_",i,j,sep=""))
					subMx=mapply(function(X,Y,Z){
						pnumber1=Z[i]+1
						pnumber2=Z[i+1]
						pnumber3=Z[j]+1
						pnumber4=Z[j+1]
						X[pnumber1:pnumber2,pnumber3:pnumber4]=Y
						X[pnumber3:pnumber4,pnumber1:pnumber2]=t(Y)
						return(X)
						},X=subMx,Y=subMxx,Z=pnumber)
					}
				}
			}
		}
	subMx=mapply('+',subMx,diagmx)
	list_evector=lapply(1:n_positive,function(i) evector[,i])

	meanVars=list()
	for(i in 1:p){
		gr_i=seq(min(unlist(vars[[i]]$Lt)),max(unlist(vars[[i]]$Lt)),length.out=grid_x[i])
		a1=brk1[i]+1
		a2=brk1[i+1]
		subevector=evector[a1:a2]
		list_evector_i=lapply(list_evector,'[',a1:a2)
		phimx0=lapply(vars[[i]]$Lt,function(Y){t(sapply(list_evector_i,function(x){approxfun(gr_i,x)(Y)}))})
		if(i==1){phimx=phimx0}else{phimx=mapply(function(X,Y){cbind(X,Y)},X=phimx,Y=phimx0)}
		if(i==1){meanX=muList[[i]](gr_i)}else{meanX=c(meanX,muList[[i]](gr_i))}
		meanVars[[i]]=muList[[i]](gr_i)
		}

	return(list(preCov=preCov,croCov=croCov,sigma=sigma,yCov=yCov,m_bwMu=m_bwMu,m_bwCov=m_bwCov,varsDemean=varsDemean,subMx=subMx,meanVars=meanVars,phimx=phimx,meanX=meanX,n_positive=n_positive,diagmx=diagmx,evalue=evalue,evector=evector))
	}

lengthVars <- function(vars_train, subset) {
	lenEach <- sapply(vars_train, function(x) {
		if (is.list(x)) {
		x[['userBwMu']]<-NULL
		x[['userBwCov']]<-NULL
			sapply(x, length)
		}else{
			stop('Cannot subset variable')
			}
		},simplify=FALSE)
  	len <- unique(unlist(lenEach))
  	if (length(len) != 1) {
		stop('Length of variables are not the same!')
		}
	return(len)
	}

sparse_pcacomp<-function(productrel_2,varsDemean,method,FVEthreshold,evalue,n_positive,isNewSub,diagmx,phimx){
	if(length(which(isNewSub==0))==0){varsDemean_train=varsDemean}else{varsDemean_train=varsDemean[-which(isNewSub==0)]}
	if(length(which(isNewSub==0))==0){productrel_train=productrel_2}else{productrel_train=productrel_2[-which(isNewSub==0)]}
	if(length(which(isNewSub==0))==0){phimx_train=phimx}else{phimx_train=phimx[-which(isNewSub==0)]}
	if(length(which(isNewSub==0))==0){diagmx_train=diagmx}else{diagmx_train=diagmx[-which(isNewSub==0)]}
	if(method=='FVE'){ratio=cumsum(evalue[which(evalue>0)])/sum(evalue[which(evalue>0)]);n_pcomp=min(which(ratio>FVEthreshold))}	
	if(method=='AIC'){
		logli=array(0,c(n_positive))
		for(i in 1:n_positive){
			if(i==1){psuVec=mapply(function(X,Y,Z){Z-as.matrix(X[1:i,])*Y[1:i]},X=phimx_train,Y=productrel_train,Z=varsDemean_train)}else{psuVec=mapply(function(X,Y,Z){Z-as.matrix(t(X[1:i,]))%*%as.matrix(Y[1:i])},X=phimx_train,Y=productrel_train,Z=varsDemean_train)}		
			logli[i]=sum(mapply(function(X,Y){t(Y)%*%diag(1/diag(X))%*%Y/2},X=diagmx_train,Y=psuVec))+i
			}
		ratio_vec=(max(logli)-logli)/diff(range(logli))
		n_pcomp=min(which(ratio_vec>0.95))
		#ratio_vec=abs(logli/min(logli))
		#n_pcomp=min(which(ratio_vec<1.2))
		#n_pcomp=which.min(logli>min(logli)*1.1)
		}

	if(method=='BIC'){
		logli=array(0,c(n_positive))
		for(i in 1:n_positive){
			if(i==1){psuVec=mapply(function(X,Y,Z){Z-as.matrix(X[1:i,])*Y[1:i]},X=phimx_train,Y=productrel_train,Z=varsDemean_train)}else{psuVec=mapply(function(X,Y,Z){Z-as.matrix(t(X[1:i,]))%*%as.matrix(Y[1:i])},X=phimx_train,Y=productrel_train,Z=varsDemean_train)}		
			logli[i]=sum(mapply(function(X,Y){t(Y)%*%diag(1/diag(X))%*%Y/2},X=diagmx_train,Y=psuVec))+i*log(length(unlist(varsDemean_train)))/2
			}
		ratio_vec=(max(logli)-logli)/diff(range(logli))
		n_pcomp=min(which(ratio_vec>0.95))
		#ratio_vec=abs(logli/min(logli))
		#n_pcomp=min(which(ratio_vec<1.2))
		#n_pcomp=which.min(logli>min(logli)*1.2)
		}
	n_pcomp
	}


dense_pcacomp<-function(method,MxAll,Mxscore,evector,evalue,isNewSub,sigmaVar,n_positive,FVEthreshold){
	if(method=='FVE'){ratio=cumsum(evalue[which(evalue>0)])/sum(evalue[which(evalue>0)]);n_pcomp=min(which(ratio>FVEthreshold))}
	if(method=='AIC'){
		if(length(which(isNewSub==0))==0){MxAll_train=MxAll}else{MxAll_train=MxAll[-which(isNewSub==0),]}
		if(length(which(isNewSub==0))==0){Mxscore_train=Mxscore}else{Mxscore_train=Mxscore[-which(isNewSub==0),]}
		logli=array(0,c(n_positive))
		for(i in 1:n_positive){
			#sigmaVar_corr=sigmaVar+sum(evalue[-c(1:i)])
			coordi=MxAll_train-Mxscore_train[,1:i]%*%t(evector[,1:i])	
			logli[i]=sum(diag(coordi%*%diag(1/sigmaVar)%*%t(coordi)))/2+i
			}
		ratio_vec=(max(logli)-logli)/diff(range(logli))
		n_pcomp=min(which(ratio_vec>0.95))
		#ratio_vec=abs(logli/min(logli))
		#n_pcomp=min(which(ratio_vec<1.3))
		}
	if(method=='BIC'){
		if(length(which(isNewSub==0))==0){MxAll_train=MxAll}else{MxAll_train=MxAll[-which(isNewSub==0),]}
		if(length(which(isNewSub==0))==0){Mxscore_train=Mxscore}else{Mxscore_train=Mxscore[-which(isNewSub==0),]}
		total_p=length(MxAll_train)
		logli=array(0,c(n_positive))
		for(i in 1:n_positive){
			coordi=MxAll_train-Mxscore_train[,1:i]%*%t(evector[,1:i])	
			logli[i]=sum(diag(coordi%*%diag(1/sigmaVar)%*%t(coordi)))/2+i*log(total_p)/2
			}
		ratio_vec=(max(logli)-logli)/diff(range(logli))
		n_pcomp=min(which(ratio_vec>0.95))
		#ratio_vec=abs(logli/min(logli))
		#n_pcomp=min(which(ratio_vec<1.3))
		}
	return(n_pcomp)
	}


dense_pcacomp2<-function(method,MxAll,Mxscore,evector,evalue,isNewSub,sigmaVar,FVEthreshold){
	if(method=='FVE'){ratio=cumsum(evalue[which(evalue>0)])/sum(evalue[which(evalue>0)]);n_pcomp=min(which(ratio>FVEthreshold))}
	if(method=='AIC'){
		MxAll_train=MxAll[-which(isNewSub==0),]
		Mxscore_train=Mxscore[-which(isNewSub==0),]
		i=1
		coordi=MxAll_train-Mxscore_train[,1:i]%*%t(evector[,1:i])	
		logli_1=sum(coordi%*%diag(1/sigmaVar)%*%t(coordi))/2+i
		logli_0=logli_1+1
		logli_diff=logli_1-logli_0
		while(logli_diff<0){
			i=i+1
			logli_0=logli_1
			coordi=MxAll_train-Mxscore_train[,1:i]%*%t(evector[,1:i])	
			logli_1=sum(coordi%*%diag(1/sigmaVar)%*%t(coordi))/2+i
			logli_diff=logli_1-logli_0
			}
		n_pcomp=i-1
		}
	if(method=='BIC'){
		MxAll_train=MxAll[-which(isNewSub==0),]
		Mxscore_train=Mxscore[-which(isNewSub==0),]
		total_p=length(MxAll_train)
		i=1
		coordi=MxAll_train-Mxscore_train[,1:i]%*%t(evector[,1:i])	
		logli_1=sum(coordi%*%diag(1/sigmaVar)%*%t(coordi))/2+i
		logli_0=logli_1+1
		logli_diff=logli_1-logli_0
		while(logli_diff<0){
			i=i+1
			logli_0=logli_1
			coordi=MxAll_train-Mxscore_train[,1:i]%*%t(evector[,1:i])	
			logli_1=sum(coordi%*%diag(1/sigmaVar)%*%t(coordi))/2+i*log(total_p)/2
			logli_diff=logli_1-logli_0
			}
		n_pcomp=i
		}
	return(n_pcomp)
	}

