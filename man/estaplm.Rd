
\name{estaplm}
\alias{estaplm}
\title{Bayesian estimation for additive paritial linear models}
\usage{
estaplm(y, X, Z = NULL, nKnot = rep(20,ncol(X)),
        mcmc.n, mcmc.burnin)
}
\arguments{
\item{y}{response vector.}
\item{X}{matrix or data frame; predictor variables for nonlinear components.}
\item{Z}{matrix or data frame; predictor variables for linear components. The default is NULL.}
\item{nKnot}{integer vector; the number of knots for smooth functions. The default is 20 for each predictor in X.}
\item{mcmc.n}{integer; the number of MCMC samples after burnin. The default is 2000.}
\item{mcmc.burnin}{integer; the number of burnin-iterations. The default is 100.}
}
\description{
This function is to estimate an additive partial linear model under a certain model formulation. An additive model can also be treated since it is a subclass of an additive partial linear model.
}
\value{
\item{grid.of.x}{grid of x for function evaluation.}
\item{function.estimate}{MCMC estimates of the smooth functions.}
\item{parameter.estimate}{MCMC estimates of the parameters other than the smooth functions.}
\item{function.summary}{posterior summary for the smooth functions.}
\item{parameter.summary}{posterior summary for the parameters other than the smooth functions.}
}
\seealso{
\code{\link{bmsaplm}}
}
\author{
Seonghyun Jeong
}
\examples{
data(betacaro)
dataset=betacaro[betacaro$logbetacaro!=-Inf,]
y=dataset$logbetacaro
X=dataset[,"age"]
Z=dataset[,c("fiber","bmi","sex","smokestat","vituse")]
# X and Z are chosen by model selection; see the function bmsaplm.

est1=estaplm(y,X,Z)

# est2=estaplm(y,X)     # without linear components.

summary(est1)

plot(est1)
}

