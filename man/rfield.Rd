\name{rfield}
\alias{rfield}
\title{
Simulate data from a 2-D field with non-radial basis
}
\description{
Currently, this function can a field simulate from a regular grid of basis functions with normal errors.
}
\usage{
rfield(gridlen, coefs, shape, rate, scale, meanlog, sdlog,
        nbasis, type = c("gamma", "weibull", "lognormal"), rev = FALSE,
        copula = FALSE, copulaType = NULL, param = NULL,
        errtype = c("gaussian"), errvar = 1)
}
\arguments{
  \item{gridlen}{Length of regular grid.}
  \item{coefs}{Numeric vector of coefficients, including intercept term (if no intercept, set first entry to 0). Length of coefs must be equal to nbasis^2.}
  \item{shape}{Shape parameter of gamma of weibull distributions.}
  \item{rate}{Rate parameter of gamma distribution.}
  \item{scale}{Scale parameter of weibull distribution.}
  \item{meanlog}{Mean of the normal distribution on the log scale}
  \item{sdlog}{Standard deviation of the normal distribution on the log scale}
  \item{nbasis}{Number of basis functions to create in 1 dimension (this gets squared if making a field).}
  \item{rev}{Logical. Should the direction of the basis functions be reverse? In 1-D, this correspond to letting y=f(-x).}
  \item{copula}{Logical. Should 2-D bases be made up of dependent 1-D components?}
  \item{copulaType}{A character vector specifying which type of copula should be used to model the dependence when copula = TRUE. Options include "clayton", "frank", "gumbel", "amh", and "joe". These are Archimedean copulas implemented in the copula package.}
  \item{param}{Numeric value of the copula dependence parameter.}
  \item{errtype}{What sort of error to add to the basis? Currently, only gaussian errors are available.}
  \item{errvar}{Numeric specifying the variance of the error distribution.}
}
\details{
Currently no support for irregular grid simulation (aside by selecting some coefficients to be 0).
}
\value{
Matrix corresponding a linear combination of basis functions.
}
\author{\email{hbk5086@psu.edu}}
\seealso{\code{\link[copula]{archmCopula}}}
\examples{
library(VGAM)
set.seed(1)

nbasis <- 40
coefs <- VGAM::rgpd(nbasis^2+1, location = 1, scale = 10, shape = 0.75)
coefs[sample(seq_along(coefs), round(length(coefs)*.3), replace = FALSE)] <- 0

b <- create.weibull.field(80, 1.5, 3, nbasis, rev = TRUE, reduce = TRUE)
pure <- rfield(80, coefs = coefs, shape = 1.5, scale = 3, nbasis = nbasis,
                type = "weibull", rev = TRUE, errvar = 0)
noisy <- rfield(80, coefs = coefs, shape = 1.5, scale = 3, nbasis = nbasis,
                type = "weibull", rev = TRUE, errvar = 1)

par(mfrow = c(1,3))
contour(b)
contour(pure)
contour(noisy)

coefs2 <- VGAM::rgpd(nbasis^2+1, location = 0, scale = 5, shape = 0.65)
b2 <- create.lognormal.field(100, .5, 3, nbasis, rev=TRUE,reduce=TRUE, copula = TRUE, copulaType = "frank", param = -4)
pure2 <- rfield(80, coefs = coefs2, meanlog = .5, sdlog = 3, nbasis = nbasis,
                type = "lognormal", rev = TRUE, copula = TRUE, copulaType = "frank", param = -4, errvar = 0)
noisy2 <- rfield(80, coefs = coefs2, meanlog = .5, sdlog = 3, nbasis = nbasis,
                type = "lognormal", rev = TRUE,  copula = TRUE, copulaType = "frank", param = -4, errvar = .5)

image(b2^.3, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE)
image(pure2^.3, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE)
image(noisy2^.3, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE)
}
