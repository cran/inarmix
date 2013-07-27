\name{GenerateMixData}
\alias{GenerateMixData}
\title{
Generates simulated data according to a specified INAR mixture model.
}
\description{
This function generates data according to a user-specified INAR
mixture model. It returns a data frame which can be used immediately
by inarmix.
}
\usage{
GenerateMixData(m, coefs, autocorr, scale=NULL, mix, design.mat,
                return.labels=F,poisson=F)
}
\arguments{
  \item{m}{
a numeric scalar (the number of subjects)
}
  \item{coefs}{
a numeric vector
}
  \item{autocorr}{
a numeric scalar
}
  \item{scale}{
a numerica scalar
}
  \item{mix}{
a numeric vector
}
  \item{design.mat}{
a numeric matrix (the design matrix)
}
  \item{return.labels}{
an indicator of whether or not to return the class labels
}
  \item{poisson}{
an indicator of whether or not the data should have 
marginal Poisson distributions.
}
}
\value{
A data frame which contains the response and a subject identifier.
The other columns contain the data entered from the design.mat argument.
}
\author{
Nicholas Henderson
}
\examples{
XX <- cbind(rep(1,9),c(0:8)/4)
colnames(XX) <- c("const","time")
coefs <- rbind(c(-.2,0),c(1.2,.3))
autocorr <- c(.2,.2)
scale <- c(1.5,1.5)
mix.prop <- c(.8,.2)

testdat <- GenerateMixData(1000,coefs,autocorr,scale,mix.prop,XX)
}
\keyword{datagen}
