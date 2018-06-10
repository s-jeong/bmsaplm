bmsaplm=function(y, X, Z=NULL, nKnot=rep(20,ncol(X)), pilot.n=300, main.n=2000, pilot.burnin=50, main.burnin=50){
	this.call=match.call()
	pilot.num=pilot.n+pilot.burnin
	main.num=main.n+main.burnin
	y=as.vector(y)
	X=as.matrix(X)
	n=length(y)
	p=ncol(X)
	p.ztgeo=rep(0.4,p)
	colnamesX=colnames(X)
	if(is.null(Z)){
		colnamesZ=NULL
		r=0
	}else{
		Z=as.matrix(Z)
		colnamesZ=colnames(Z)
		if(is.null(colnamesZ)) colnamesZ=paste("Z",1:r,sep="")
		Z=sweep(Z,2,colMeans(Z))
		r=ncol(Z)
	}
	if(is.null(colnamesX)) colnamesX=paste("X",1:p,sep="")
	if(r!=0){
		if(!all(c(nrow(X),nrow(Z))==n)) stop("The number of observations in y, X, and Z should be the same.")
	} else {
		if(!(nrow(X)==n)) stop("The number of observations in y and X should be the same.")
	}
	if(!all(nKnot>=1)) stop("Use at least one knot for each covariate")
	if(main.n<2 | main.burnin<1 | pilot.n<2 | pilot.burnin<1) stop("Use appropriate values for the number of MCMC samples.")
	numkn=nKnot
	makelist.col=function(x) if(!is.null(ncol(x))){lapply(seq_len(ncol(x)),function(i) x[,i])}else{list(x)}
	list.X=makelist.col(X)
	knM=mapply(function(x,y) quantile(unique(x),ppoints(y,a=0)),list.X,numkn,SIMPLIFY=F)
	uncent.WstarM=CombBasis(X,knM)
	WstarM=sweep(uncent.WstarM,2,colMeans(uncent.WstarM))
	WstarM.Zi=cbind(WstarM,Z)
	nfix=sum(numkn)+p+r
	cor.vec=c()
	for(j in 1:p){
		temp.W=WstarM.Zi[,(sum(numkn[0:(j-1)])+j):(sum(numkn[1:j])+j)]
		corZ=cor(temp.W)
		cor.vec[j]=max(corZ[lower.tri(corZ)])
	}
	if(any(cor.vec>1-1e-7)){
		stop(paste("Use the fewer number of knot candidates for ",
			paste(colnamesX[c(cor.vec>1-1e-7)],collapse=", "),
			". The current valeus of nKot are ",paste(numkn,collapse=", "),".",sep=""))
	}
	delta=lapply(numkn,function(x)rep(0,x+1))
	for(dd in 1:p) {delta[[dd]][1]=1;delta[[dd]][sample(2:length(delta[[dd]]),1)]=1}
	delta.var=c(rep(2,p),rep(1,r))

	temp.list.delta=list(list(0,1),
				list(c(0,0),c(1,0),c(0,1),c(1,1)),
				list(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,0),c(0,0,1),c(1,0,1),c(0,1,1),c(1,1,1)),
				list(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),c(1,1,0,0),c(0,0,1,0),c(1,0,1,0),c(0,1,1,0),c(1,1,1,0),
					c(0,0,0,1),c(1,0,0,1),c(0,1,0,1),c(1,1,0,1),c(0,0,1,1),c(1,0,1,1),c(0,1,1,1),c(1,1,1,1)))
	Fset.delta.comb=matrix(,length(unlist(delta)),pilot.num)
	v.delta=as.vector(Eta(delta,delta.var,p,r))
	log.BF.cur=LogBF(v.delta,WstarM.Zi,y,n,-3/4)
	cat("Pilot chain:",'\n')
	cat("0% =================== 50% =================== 100%",'\n')
	for(iter in 1:pilot.num){
		delta.iter=PilotSampleDeltaGlobalUpdateVec(WstarM.Zi,y,delta,delta.var,temp.list.delta,numkn,n,p,r,log.BF.cur,p.ztgeo)
		delta=delta.iter[[1]]
		log.BF.cur=delta.iter[[2]]
		Fset.delta.comb[,iter]=unlist(delta)
		if(iter%%round(max(1,pilot.num/50))==0) cat(paste(rep("+",round(50*iter/pilot.num)),collapse=""),"\r")
	}
	cat("",'\n');cat("",'\n')
	Fset.delta=lapply(split(as.data.frame(Fset.delta.comb),rep(1:p,mapply(length,delta))),function(x)t(as.matrix(x)))
	set.delta=lapply(Fset.delta,function(x) x[-(1:pilot.burnin),])
	sum.delta=sapply(set.delta,function(x)apply(x[,-1],1,sum))
	mean.delta0=sapply(set.delta,function(x)mean(x[,1]))
	mean.delta0=(mean.delta0-0.5)*0.99+0.5
	temp.m.delta=apply(sum.delta,2,mean)
	temp.var.delta=apply(sum.delta,2,function(x)max(0.01,var(x)))
	m.delta=var.delta=c()
	for(j in 1:p){
		findpar=function(x) {
			y=c()
			mean.d=sum((1:numkn[j])*dnorm(1:numkn[j],x[1],abs(x[2]))/sum(dnorm(1:numkn[j],x[1],abs(x[2]))))
			var.d=sum((1:numkn[j])^2*dnorm(1:numkn[j],x[1],abs(x[2]))/sum(dnorm(1:numkn[j],x[1],abs(x[2]))))-mean.d^2
			y[1]=mean.d-temp.m.delta[j]; y[2]=var.d-temp.var.delta[j]
			return(y)
		}
		sol=nleqslv(c(1,1),findpar)$x
		m.delta[j]=sol[1]
		var.delta[j]=sol[2]^2
	}
	delta.var=rep(1,p+r)
	Fset.sigma.sq=Fset.qv=Fset.int=c()
	Fset.delta.comb=matrix(,length(unlist(delta)),main.num)
	Fset.delta.var=matrix(NA,main.num,p+r)
	v.delta=Eta(delta,delta.var,p,r)
	log.BF.cur=LogBF(v.delta,WstarM.Zi,y,n,-3/4)
	cat("Main chain:",'\n')
	cat("0% =================== 50% =================== 100%",'\n')
	for(iter in 1:main.num){
		MCMC.one=MCMCOneIteration(WstarM.Zi,y,delta,delta.var,temp.list.delta,numkn,n,p,r,log.BF.cur,p.ztgeo,m.delta,var.delta,mean.delta0,-3/4)
		delta=MCMC.one[[1]]
		delta.var=MCMC.one[[2]]
		log.BF.cur=MCMC.one[[3]]
		Fset.delta.comb[,iter]=unlist(delta)
		Fset.delta.var[iter,]=delta.var
		if(iter%%round(max(1,main.num/50))==0) cat(paste(rep("+",round(50*iter/main.num)),collapse=""),"\r")
	}
	cat("",'\n');cat("",'\n')
	Fset.delta=lapply(split(as.data.frame(Fset.delta.comb),rep(1:p,mapply(length,delta))),function(x)t(as.matrix(x)))
	set.delta=lapply(Fset.delta,function(x) x[-(1:main.burnin),])
	set.delta.var=as.matrix(Fset.delta.var[-(1:main.burnin),])
	p.delta=lapply(set.delta,function(x) apply(as.matrix(x),2,mean))
	p.delta.var=matrix(as.vector(sapply(0:2,function(k)apply(set.delta.var==k,2,mean))),nrow=3,byrow=T)
	colnames(p.delta.var)=c(colnamesX,colnamesZ)
	rownames(p.delta.var)=c("Zero","Linear","Nonlinear")
	gamma.est=apply(p.delta.var,2,which.max)-1
	structure(list(inclusion.probability=p.delta.var,latent.variable.estimate=gamma.est,call=this.call),class="bmsaplm")
}

