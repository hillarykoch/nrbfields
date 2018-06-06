\name{create}
\alias{create.weibull.basis}
\alias{create.gamma.basis}
\alias{create.weibull.field}
\alias{create.gamma.field}
\title{
Create 1- or 2-D non-radial basis
}
\description{
This collection of functions can be used to generate fields of asymmetric basis functions. Basis functions can be reversed in direction, depending on the data being examined.
}
\usage{
create.weibull.basis(gridlen, shape, scale, nbasis, rev = FALSE)
create.gamma.basis(gridlen, shape, rate, nbasis, rev = FALSE)
create.weibull.field(gridlen, shape, scale, nbasis, rev=FALSE, reduce = FALSE)
create.gamma.field(gridlen, shape, rate, nbasis, rev = FALSE, reduce = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gridlen}{Length of regular grid.}
  \item{shape}{Shape parameter of gamma of weibull distributions.}
  \item{scale}{Scale parameter of weibull distribution.}
  \item{rate}{Rate parameter of gamma distribution.}
  \item{nbasis}{Number of basis functions to create in 1 dimension (this gets squared if making a field).}
  \item{rev}{Logical. Should the direction of the basis functions be reverse? In 1-D, this correspond to letting y=f(-x).}
  \item{reduce}{Logical. Should list of 2-D bases be summed into one field?}
}
\details{
Currently, no support for generating irregular grids.
}
\value{
For create.basis, a basis gridlen by nbasis basis matrix is returned. For create field, if reduce = TRUE, a nbasis by nbasis matrix is returned. Otherwise, a list of length nbasis^2, each element containing a gridlen by gridlen matrix corresponding to one basis function, is returned.
}
\author{\email{hbk5086@psu.edu}}
\examples{
2dfield <- create.gamma.field(100, 3, .5, 10, reduce = TRUE)
contour(2dfield)
}