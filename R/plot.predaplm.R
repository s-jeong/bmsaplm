plot.predaplm=function(x,kde=FALSE,...){
	if(kde==TRUE){
		plot(density(x$posterior.predictive.values),xlab="",ylab="",main="")
	} else {
		hist(x$posterior.predictive.values,prob=T,main="",xlab="",ylab=""); box()
	}
	mtext("Density",side=2,line=2.3,cex=0.9)
	mtext("Posterior predictive values",side=1,line=2.3,cex=0.9)
}
