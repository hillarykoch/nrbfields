\name{call_corners}
\alias{call_corners}
\title{
Quickly find corners in contact matrix
}
\description{
Use the Canny edge detector to find edges in a smoothed contact matrix, and determine corners based on the results.
}
\usage{
call_corners(x, alpha, pow=0.2, t=3, t1=NULL, t2=NULL, sigma=2)
}
\arguments{
  \item{x}{Numeric square contact matrix}
  \item{alpha}{threshold adjustment factor for the Canny edge detector}
  \item{pow}{Power to raise elements of the contact matrix to for image construction. Defaults to 1.}
  \item{t}{Diffusion parameter for the initial contact matrix. Defaults to 1.}
  \item{t1, t2}{Threshold for weak and strong edges.}
  \item{sigma}{Degree of smoothing for the image gradient. Defaults to 2.}
}
\details{
This function calls \link[imager]{cannyEdges} from the imager package, and arguments alpha, t1, t2, and sigma are exactly those described in the documentation there.
}
\value{A binary, upper triangular matrix that is the same size as x. Corner loctions are indicated by 1.}
\author{\email{hbk5086@psu.edu}}
\examples{
set.seed(1)
nbasis <- 15

b <- create.lognormal.field(100, .5, 3, nbasis, rev = TRUE, reduce = FALSE,
                            copula = TRUE, copulaType = "frank", param = -4)
b <- b[-sample(1:nbasis^2, 75, replace = FALSE)]

coefs <- rexp(nbasis^2-50) + 1
fld <- lapply(seq_along(b), function(X) b[[X]]*coefs[X])
fld <- Reduce('+', fld) + matrix(rnorm(100^2, 0, .008), nrow = 100)

fld <- fld + abs(min(fld))
fld[lower.tri(fld)] <- 0


corners <- call_corners(fld, alpha = .4, pow = 0.2)

par(mfrow = c(1,2))
image(fld^.4, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE, main = "image")
image(corners, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE, main = "called corners")

}