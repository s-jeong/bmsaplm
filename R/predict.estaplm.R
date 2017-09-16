predict.estaplm=function(object,x,z=NULL,...){
	this.call=match.call()
	if(!all(c(is.vector(x),any(c(is.null(z),is.vector(z)))))) stop("x and z should be vectors (or constants).")
	if(!all(c(length(object$function.summary)==length(x),(nrow(object$parameter.summary)-2)==length(z)))) stop("The length of x or z is wrong.")
	fn.eval=mapply(function(w,t,s) apply(w,2,function(m) approx(t,m,s)$y ), object$function.estimate,object$grid.of.x,as.list(x))
	mu=as.matrix(object$parameter.estimate[,-ncol(object$parameter.estimate)])%*%as.matrix(c(1,z))+rowSums(fn.eval)
	y.tilde=rnorm(length(mu),mu,sqrt(object$parameter.estimate[,ncol(object$parameter.estimate)]))
	summary=c(quantile(y.tilde,c(0.025,0.5,0.975)),mean(y.tilde),var(y.tilde))
	names(summary)=c("2.5% percentile","Median","97.5% percentile","Mean","Variance")
	structure(list(posterior.predictive.values=y.tilde, posterior.predictive.summary=summary,call=this.call),class="predaplm")
}
