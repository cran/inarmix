\name{diagnose}
\alias{diagnose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostics for the model fit.
}
\description{
Provides diagnostics of the mixture model fit.
}
\usage{
diagnose(results)
}
\arguments{
\item{results}{
An object of class "inarmix"
} 
}
\value{
An object of class \code{"diagnose.inarmix"} which is a list containing the 
following components
\item{em.converged}{
An indication of whether or not the "EM algorithm" converged.
}
\item{niter}{
The number of iterations for required for the EM algorithm to converge.
}
\item{nclasses}{
The number of latent classes used in the model fit.
}
\item{loglikhood}{
The values of the log-likelihood function at each EM iteration.
}
\item{ConvMat}{
This is a matrix containing several pieces of information.
The columns with GEE in the title indicate if the
algorithms for updating the parameters within an EM
iteration have converged. The column with l[i] just
contains the log-likelihood value. The column with
the title "||Psi^2||" gives the norm of the "global
estimating equation" - this should be near zero if the
algorithm converged.
}
}
\author{
Nicholas Henderson
}
\examples{
XX <- cbind(rep(1,9),c(0:8)/4)
colnames(XX) <- c("const","time")
coefs <- rbind(c(-.2,0),c(1.2,.3))
alpha <- c(.2,.2)
scale <- c(2,2)
mix.prop <- c(.8,.2)

testdat <- GenerateMixData(200,coefs,alpha,scale,mix.prop,XX)
testfit <- inarmix(y~time,nclasses=2,id=subject,data=testdat,maxiter=3)

diagnose(testfit)
}
\keyword{print}
