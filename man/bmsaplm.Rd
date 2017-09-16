
\name{bmsaplm}
\alias{bmsaplm}
\title{Bayesian model selection for additive paritial linear models}
\usage{
bmsaplm(y, X, Z = NULL, nKnot = rep(20,ncol(X)),
        pilot.n, main.n, pilot.burnin, main.burnin)
}
\arguments{
\item{y}{response vector}
\item{X}{matrix or data frame; predictor variables for nonlinear components.}
\item{Z}{matrix or data frame; predictor variables for linear components. The default is NULL.}
\item{nKnot}{integer vector; the number of knots for smooth functions. The default is 20 for each predictor in X.}
\item{pilot.n}{integer; the number of MCMC samples for a pilot chain after burnin. The default is 300.}
\item{main.n}{integer; the number of MCMC samples for a main chain after burnin. The default is 2000.}
\item{pilot.burnin}{integer; the number of burnin-iterations for a pilot chain. The default is 50.}
\item{main.burnin}{integer; the number of burnin-iterations for a main chain. The default is 50.}
}
\description{
This function is to perform model selection for an additive partial linear model using spline basis search. An additive model can also be treated since it is a subclass of an additive partial linear model.
}
\value{
\item{inclusion.probability}{inclusion probability of each predictor.}
\item{latent.variable.estimate}{estimates of the latent variables from the inclusion probabilities.}
}
\seealso{
\code{\link{estaplm}}
}
\author{
Seonghyun Jeong
}
\examples{
data(betacaro)
dataset=betacaro[betacaro$logbetacaro!=-Inf,]
y=dataset$logbetacaro
X=dataset[,c("age","bmi","calories","fat","fiber","alcohol","chol","betadiet","retdiet")]
Z=dataset[,c("sex","smokestat","vituse")]

# bmsaplm(y,X,Z)   # this does not work.

bms1=bmsaplm(y,X,Z,nKnot=c(20,20,20,20,20,8,20,20,20))   # this works.

# bms2=bmsaplm(y,X,nKnot=c(20,20,20,20,20,8,20,20,20))     # without linear components.

summary(bms1)
}
