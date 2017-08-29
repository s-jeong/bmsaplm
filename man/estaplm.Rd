
\name{estaplm}
\alias{estaplm}
\title{Bayesian estimation for additive paritial linear models}
\usage{
estaplm=function(y,X,Z=NULL,nKnot=rep(20,ncol(X)),mcmc.n=2000,mcmc.burnin=100)
}
\arguments{
\item{y}{a vector; a response variable.}
\item{X}{a matrix or a data frame; predictor variables for nonlinear components.}
\item{Z}{a matrix or a data frame; predictor variables for linear components; the default is NULL.}
\item{nKnot}{a vector; the number of knots for X; the default is 20 for each of X.}
\item{mcmc.n}{an integer; the number of samples after burnin; the default is 2000.}
\item{mcmc.burnin}{an integer; the number of burnin-periods; the default is 100.}
}
\description{
Find the logarithm of a Bayes factor
}
\examples{
data(betaplasma)
dataset=data.frame(data.matrix(betaplasma))
dataset=data.frame(log(dataset$betacaro),dataset)
dataset=dataset[dataset$log.dataset.betacaro.!=-Inf,]
y=dataset$log.dataset.betacaro.
X=dataset[,"age"]
Z=dataset[,c("fiber","bmi","sex","smokestat","vituse")]
# X and Z are chosen by model selection; see the function bmsaplm.

estaplm(y,X,Z)

estaplm(y,X)     # only with the nonlinear component, Age.
}
