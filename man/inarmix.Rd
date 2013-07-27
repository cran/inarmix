\name{inarmix}
\alias{inarmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finite mixture model for longitudinal count data.
}
\description{
Estimates model parameters of a finite mixture model.
Appropriate for longitudinal data when the responses are counts
and when the correlation structure is assumed to be AR(1).
}
\usage{
inarmix(formula,nclasses=1,id,data,initvals=NULL,maxiter=200,stoptol=1e-7,
        num.restarts=1,time=NULL,dropthresh=.01)
}
\arguments{
\item{formula}{a formula expression used to specify the regression model}
  \item{nclasses}{
The number of components in the finite mixture model.
}
  \item{id}{
the name of the variable which identifies the individual subjects.
}
  \item{data}{
a data frame
}
  \item{initvals}{
initial estimates of the parameters (optional). This should be a list
of the form list(coef=,autocorr=,scale=,mix.prop=)
}
  \item{maxiter}{
The maximum number of EM iterations to be performed (optional).
}
\item{stoptol}{
tolerance level which determines convergence. The default is 1e-7.
}
\item{num.restarts}{
The number of runs. Each run has a random starting value for the 
parameters.
}
\item{time}{
the name of the variable which indicates time. When left blank, the 
data are assumed to have the correct time-ordering. 
}
\item{dropthresh}{
The threshold at which one class is dropped from the model.
If the estimated proportion for a class drops below this level,
it is removed from the fitting procedure.
}
}
\value{
An object of class \code{"inarmix"} which is a list containing at least
the following components
\item{coefficients}{A matrix of estimated regression coefficients. Each row 
contains the coefficients for one class.}
\item{mix.prop}{The estimated class-membership probabilities.}
\item{post.probs}{A nclasses x num.subjects matrix. The posterior probabilities
of class-membership for each subject and class.}
\item{loglikhood}{The final value of the log-likelihood.}
\item{niter}{The number of iterations required for convergence.}
\item{cov.mat}{The variance covariance matrix of the parameter estimates.}
\item{call}{the matched call}
\item{nclasses}{
The number of classes in the final model fit.
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Nicholas Henderson
}
\examples{
set.seed(4297)

############################################################
#### Simulate data from a two class model

XX <- cbind(rep(1,9),c(0:8)/4)
colnames(XX) <- c("const","time")
beta <- rbind(c(-.2,0),c(1.2,.3))
### this means that for group 1: (beta_{0},beta_{1}) = (-.2,0)
### and for group 2: (beta_{0},beta_{1}) = (1.2,.3)
autocorr <- c(.2,.2)
scale <- c(2,2)
mix.prop <- c(.8,.2)  ## proportion in group 1 is .8

testdat <- GenerateMixData(500,beta,autocorr,scale,mix.prop,XX)
testdat[1:5,]


########################################################################
####  Fit a linear curve with two classes (with a maximum of 4 iterations)

twoclassfit <- inarmix(y~time,nclasses=2,id=subject,data=testdat,maxiter=4)
summary(twoclassfit)

diagnose(twoclassfit)

#############################################################
#####  Fit the same model with specified starting values.

inpars <- list()
inpars$coef <- rbind(c(-.5,.1),c(.5,0))
inpars$autocorr <- rep(.3,2)
inpars$scale <- rep(2,2)
inpars$mix.prop <- c(.6,.4)

twoclassfit2 <- inarmix(y~time,nclasses=2,id=subject,data=testdat,initvals=inpars,
                        maxiter=4)
summary(twoclassfit2)


###############################################################
### Try fitting a one class model with the same data
oneclassfit <- inarmix(y~time,nclasses=1,id=subject,data=testdat)
summary(oneclassfit)


#########################################################################
####  Fit a two class model with multiple starts
\dontrun{
testfit_multi <- inarmix(y~time,nclasses=2,id=subject,data=testdat,num.restarts=3)
summary(testfit_multi)

#### Look at final log-likelihood values for each restart
testfit_multi$reploglik

#### Look at final parameter estimates for each restart
testfit_multi$finalvals

#############################################################################
###########  Simulate data from a four class model

XX <- cbind(rep(1,9),seq(0,2,by=.25))
colnames(XX) <- c("const","time")
beta <- rbind(c(-.4,-.1),c(1.4,-.6),c(0,.7),c(1.4,0))
autocorr <- rep(.2,4)
scale <- rep(1.5,4)
mix.prop <- c(.5,.25,.15,.1)

testdat4 <- GenerateMixData(1000,beta,autocorr,scale,mix.prop,XX)

### Fit a four class model

testfit_four <- inarmix(y~time,nclasses=4,id=subject,data=testdat4,maxiter=5)
summary(testfit_four)
}

}
\keyword{models}
\keyword{regression}
