
\name{bmsaplm}
\alias{bmsaplm}
\title{Bayesian model selection for additive paritial linear models}
\usage{
bmsaplm(y,X,Z=NULL,nKnot=rep(20,ncol(X)),pilot.n=300,main.n=2000,pilot.burnin=50,main.burnin=50)
}
\arguments{
\item{y}{a vector; a response variable.}
\item{X}{a matrix or a data frame; predictor variables for nonlinear components.}
\item{Z}{a matrix or a data frame; predictor variables for linear components; the default is NULL.}
\item{nKnot}{a vector; the number of knots for X; the default is 20 for each of X.}
\item{pilot.n}{an integer; the number of samples for a pilot chain after burnin; the default is 300.}
\item{main.n}{an integer; the number of samples for a main chain after burnin; the default is 2000.}
\item{pilot.burnin}{an integer; the number of burnin-periods for a pilot chain; the default is 50.}
\item{main.burnin}{an integer; the number of burnin-periods for a main chain; the default is 50.}
}
\description{
Bayesian model selection for additive paritial linear models
}
\examples{
data(betaplasma)
dataset=data.frame(data.matrix(betaplasma))
dataset=data.frame(log(dataset$betacaro),dataset)
dataset=dataset[dataset$log.dataset.betacaro.!=-Inf,]
y=dataset$log.dataset.betacaro.
X=dataset[,c("age","bmi","calories","fat","fiber","alcohol","chol","betadiet","retdiet")]
Z=dataset[,c("sex","smokestat","vituse")]

bmsaplm(y,X,Z)   # this does not work.

bmsaplm(y,X,Z,nKnot=c(20,20,20,20,20,8,20,20,20))   # this works.

bmsaplm(y,X,nKnot=c(20,20,20,20,20,8,20,20,20))     # only with nonlinear components.
}
