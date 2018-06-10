estaplm=function(y, X, Z=NULL, nKnot=rep(20,ncol(X)), mcmc.n=2000, mcmc.burnin=100, hyper.parameter=0.4){
	this.call=match.call()
	num=mcmc.n+mcmc.burnin
	y=as.vector(y)
	X=as.matrix(X)
	n=length(y)
	p=ncol(X)
	if(any(hyper.parameter<=0 | hyper.parameter>=1)) stop("Provide a value (or a vector) that is strictly between 0 and 1 for the hyperparameter.")
	if(length(hyper.parameter)!=1 & length(hyper.parameter)!=p) stop("Provide a single value or a p-dimensional vector for the hyperparameter.")
	if(length(hyper.parameter)==1){
  	p.ztgeo=rep(hyper.parameter,p)
	} else {
	  p.ztgeo=hyper.parameter
	}
	colnamesX=colnames(X)
	if(is.null(Z)){
		colnamesZ=NULL
		r=0
	}else{
		Z=as.matrix(Z)
		colnamesZ=colnames(Z)
		colmeanZ=colMeans(Z)
		if(is.null(colnamesZ)) colnamesZ=paste("Z",1:r,sep="")
		Z=sweep(Z,2,colmeanZ)
		r=ncol(Z)
	}
	if(is.null(colnamesX)) colnamesX=paste("X",1:p,sep="")
	if(r!=0){
		if(!all(c(nrow(X),nrow(Z))==n)) stop("The number of observations in y, X, and Z should be the same.")
	} else {
		if(!(nrow(X)==n)) stop("The number of observations in y and X should be the same.")
	}
	if(!all(nKnot>=1)) stop("Use at least one knot for each covariate")
	if(mcmc.n<2 | mcmc.burnin<1) stop("Use appropriate values for the number of MCMC samples.")
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
	Fset.sigma.sq=Fset.qv=Fset.int=c()
	Fset.fixeff=matrix(0,num,nfix)
	Fset.delta.comb=matrix(,length(unlist(delta)),num)
	v.delta=Eta(delta,delta.var,p,r)
	log.BF.cur=LogBF(v.delta,WstarM.Zi,y,n,-3/4)
	cat("Estimation:",'\n')
	cat("0% =================== 50% =================== 100%",'\n')
	for(iter in 1:num){
		MCMC.est=MCMCEstIteration(WstarM.Zi,y,delta,delta.var,temp.list.delta,numkn,n,p,r,log.BF.cur,-3/4, p.ztgeo)
		delta=MCMC.est[[1]]
		log.BF.cur=MCMC.est[[3]]
		Fset.int[iter]=MCMC.est[[7]]
		Fset.qv[iter]=MCMC.est[[4]]
		Fset.sigma.sq[iter]=MCMC.est[[5]]
		Fset.fixeff[iter,MCMC.est[[2]]==1]=MCMC.est[[6]]
		Fset.delta.comb[,iter]=unlist(delta)
		if(iter%%round(max(1,num/50))==0) cat(paste(rep("+",round(50*iter/num)),collapse=""),"\r")
	}
	cat("",'\n');cat("",'\n')
	Fset.int=Fset.int-Fset.fixeff[,(ncol(Fset.fixeff)-r+1):ncol(Fset.fixeff)]%*%colmeanZ
	Fset.delta=lapply(split(as.data.frame(Fset.delta.comb),rep(1:p,mapply(length,delta))),function(x)t(as.matrix(x)))
	set.int=Fset.int[-(1:mcmc.burnin)]
	set.qv=Fset.qv[-(1:mcmc.burnin)]
	set.sigma.sq=Fset.sigma.sq[-(1:mcmc.burnin)]
	set.fixeff=Fset.fixeff[-(1:mcmc.burnin),]
	set.delta=lapply(Fset.delta,function(x) x[-(1:mcmc.burnin),])
	if(r!=0) set.beta=set.fixeff[,(nfix-r+1):nfix] else set.beta=NULL
	p.delta=lapply(set.delta,function(x) apply(as.matrix(x),2,mean))
	p.delta.f=lapply(set.delta,function(x)mean(apply(as.matrix(x),1,sum)!=0))
	output=cbind(set.int,set.beta,set.sigma.sq); colnames(output)=c("Intercept",colnamesZ,"Sigma^2")
	output.summary=round(t(apply(output,2,function(x)c(quantile(x,c(0.025,0.5,0.975)),mean(x),var(x)))),4)
	colnames(output.summary)=c("2.5% percentile","Median","97.5% percentile","Mean","Variance")
	xy = function(y1,y2,x){
		new.x  <- rep(seq(min(x)-diff(x)[1]/2,max(x)+diff(x)[1]/2,diff(x)[1]),each=2)
		new.y1 <- c(0,rep(y1,each=2),0)
		new.y2 <- c(0,rep(y2,each=2),0)
		x.rev  <- sort(new.x,decreasing=T)
		y2.tmp <- new.y2[length(new.y2):1]
		return(cbind(c(new.x,x.rev),c(new.y1,y2.tmp)))
	}
	grd=list(); for(i in 1:p) grd[[i]]=seq(min(X[,i]),max(X[,i]),length.out=500)
	fixind=c(0,cumsum(sapply(delta,length)))
	vcm.list.c=list()
	for(i in 1:p){
		uncent.WstarM=CombBasis(as.matrix(grd[[i]]),knM[i])
		WstarM=sweep(uncent.WstarM,2,colMeans(uncent.WstarM))
		vcm.list.c[[i]]=WstarM%*%t(set.fixeff[,(fixind[i]+1):fixind[i+1]])
	}
	list.fn.summary=mapply(function(y,z){
		fn.summary=t(apply(y,1,function(x) c(quantile(x,c(0.025,0.5,0.975)),mean(x),var(x))))
		colnames(fn.summary)=c("2.5% percentile","Median","97.5% percentile","Mean","Variance")
		rownames(fn.summary)=round(z,4)
		return(fn.summary)
	},vcm.list.c,grd,SIMPLIFY=F)
	names(list.fn.summary)=names(grd)=colnamesX
	structure(list(grid.of.x=grd, function.estimate=vcm.list.c, parameter.estimate=output, function.summary=list.fn.summary, parameter.summary=output.summary,call=this.call),class="estaplm")
}

