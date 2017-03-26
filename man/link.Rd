\name{link}
\alias{link}
\title{Link Functions for Aster Models}
\usage{
link(xi, fam, deriv = 0, delta)
}
\description{
Calculate link function and up to one derivative for
families known to the package.
}
\arguments{
  \item{xi}{mean value parameter value, a numeric vector.}
  \item{fam}{an object of class \code{"astfam"} produced by one of
    the family functions (see \code{\link{families}}) specifying the
    exponential family.}
  \item{deriv}{the number of derivatives wanted.  Must be either
    zero or one.}
  \item{delta}{direction in which limit is taken.  Link
    function is for family that is limit of family specified, limit
    being for distributions with canonical parameter
    \eqn{\theta + s \delta}{theta + s delta},
    the limit being as \eqn{s \to \infty}{s to infinity}.  May be missing,
    in which case \eqn{\delta = 0}{delta = 0} is assumed, which is the same
    as no limit being taken.}
}
\value{
  a list containing some of the following components:
  \item{zeroth}{the value of the link function at \eqn{\xi}{xi}, a vector
    of dimension \eqn{d}, where \eqn{d}} is the dimension of \eqn{\xi}{xi}.
  \item{first}{the value of the first derivative at \eqn{\xi}{xi}, a
    \eqn{d \times d}{d by d} matrix, where \eqn{d} is the dimension
    of \eqn{\xi}{xi} or a scalar if \eqn{\xi}{xi} is scalar.}
}
\section{Note}{
Not intended for use by ordinary users.  Provides R interface for testing
to C code called by many other functions in the package.
}
\seealso{
\code{\link{families}} and \code{\link{cumulant}}
}
\examples{
link(0.3, fam.bernoulli(), deriv = 1)
link(0.3, fam.bernoulli(), deriv = 1, delta = 1)
}
\keyword{misc}
