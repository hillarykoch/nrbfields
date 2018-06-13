\name{create}
\alias{create.weibull.basis}
\alias{create.gamma.basis}
\alias{create.lognormal.basis}
\alias{create.weibull.field}
\alias{create.gamma.field}
\alias{create.lognormal.field}
\title{
Create 1- or 2-D non-radial basis
}
\description{
This collection of functions can be used to generate fields of asymmetric basis functions. Basis functions can be reversed in direction, depending on the data being examined.
}
\usage{
create.weibull.basis(gridlen, shape, scale, nbasis, rev = FALSE, copula = FALSE, copulaType = NULL, param = NULL)
create.gamma.basis(gridlen, shape, rate, nbasis, rev = FALSE,  copula = FALSE, copulaType = NULL, param = NULL))
create.gamma.basis(gridlen, meanlog, sdlog, nbasis, rev = FALSE,  copula = FALSE, copulaType = NULL, param = NULL))
create.weibull.field(gridlen, shape, scale, nbasis, rev=FALSE, reduce = FALSE, copula = FALSE, copulaType = NULL, param = NULL))
create.gamma.field(gridlen, shape, rate, nbasis, rev = FALSE, reduce = FALSE, copula = FALSE, copulaType = NULL, param = NULL))
create.gamma.field(gridlen, meanlog, sdlog, nbasis, rev = FALSE, reduce = FALSE, copula = FALSE, copulaType = NULL, param = NULL))
}
\arguments{
  \item{gridlen}{Length of regular grid.}
  \item{shape}{Shape parameter of gamma of weibull distributions.}
  \item{scale}{Scale parameter of weibull distribution.}
  \item{rate}{Rate parameter of gamma distribution.}
  \item{nbasis}{Number of basis functions to create in 1 dimension (this gets squared if making a field).}
  \item{rev}{Logical. Should the direction of the basis functions be reverse? In 1-D, this correspond to letting y=f(-x).}
  \item{reduce}{Logical. Should list of 2-D bases be summed into one field?}
  \item{copula}{Logical. Should 2-D bases be made up of dependent 1-D components?}
  \item{copulaType}{A character vector specifying which type of copula should be used to model the dependence when copula = TRUE. Options include "clayton", "frank", "gumbel", "amh", and "joe". These are Archimedean copulas implemented in the copula package.}
  \item{param}{Numeric value of the copula dependence parameter.}
}
\details{
Currently, no support for generating irregular grids.
}
\value{
For create.basis when copula=FALSE, a gridlen by nbasis basis matrix is returned. When copula=TRUE, a gridlen by gridlen 2-D basis with only one basis function is returned (this is mostly used for being called by create.field when copula=TRUE). For create field, if reduce = TRUE, a nbasis by nbasis matrix is returned. Otherwise, a list of length nbasis^2, each element containing a gridlen by gridlen matrix corresponding to one basis function, is returned.
}
\author{\email{hbk5086@psu.edu}}
\seealso{\code{\link[copula]{archmCopula}}}
\examples{
2dfield <- create.gamma.field(100, 3, .5, 10, reduce = TRUE)
contour(2dfield)
}