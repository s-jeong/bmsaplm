
plot.estaplm=function(x,kde=FALSE,...){
	op = par(ask=TRUE)
	xy = function(y1,y2,z){
		new.x  <- rep(seq(min(z)-diff(z)[1]/2,max(z)+diff(z)[1]/2,diff(z)[1]),each=2)
		new.y1 <- c(0,rep(y1,each=2),0)
		new.y2 <- c(0,rep(y2,each=2),0)
		x.rev  <- sort(new.x,decreasing=T)
		y2.tmp <- new.y2[length(new.y2):1]
		return(cbind(c(new.x,x.rev),c(new.y1,y2.tmp)))
	}
	p=length(x$function.summary)
	colnamesX=names(x$function.summary)
	grd=x$grid.of.x
	list.fn.summary=x$function.summary
	for(i in 1:p){
		CI=xy(list.fn.summary[[i]][,1],list.fn.summary[[i]][,3],grd[[i]])
		plot(NULL,type="l",ylim=c(min(list.fn.summary[[i]][,1:3]),max(list.fn.summary[[i]][,1:3])),xlim=range(grd[[i]]),ylab="",xlab="")
		polygon(CI[-c(1,nrow(CI)/2,nrow(CI)/2+1,nrow(CI)),],col=gray(0:9/9)[8],border=F)
		lines(grd[[i]],list.fn.summary[[i]][,2],lty=1)
		mtext(paste("f(",colnamesX[i],")",sep=""),side=2,line=2.3,cex=0.9)
		mtext(colnamesX[i],side=1,line=2.3,cex=0.9)
	}
	for(i in 1:ncol(x$parameter.estimate)){
		if(kde==TRUE){
			plot(density(x$parameter.estimate[,i]),xlab="",ylab="",main="")
		} else {
			hist(x$parameter.estimate[,i],xlab="",ylab="",main="",prob=T); box()
		}
		mtext("Density",side=2,line=2.3,cex=0.9)
		mtext(colnames(x$parameter.estimate)[i],side=1,line=2.3,cex=0.9)
	}
	par(op)
}

