\name{cumulant}
\alias{cumulant}
\title{Cumulant Functions for Aster Models}
\usage{
cumulant(theta, fam, deriv = 0, delta)
}
\description{
Calculate cumulant function and up to three derivatives for
families known to the package.
}
\arguments{
  \item{theta}{canonical parameter value.}
  \item{fam}{an object of class \code{"asterfam"} produced by one of
    the family functions (see \code{\link{families}}) specifying the
    exponential family.}
  \item{deriv}{the number of derivatives wanted.  Must be nonnegative
    integer less than or equal to three.}
  \item{delta}{direction in which limit is taken.  Cumulant
    function is for family that is limit of family specified, limit
    being for distributions with parameter
    \eqn{\theta + s \delta}{theta + s delta},
    the limit being as \eqn{s \to \infty}{s to infinity}.  May be missing,
    in which case \eqn{\delta = 0}{delta = 0} is assumed, which is the same
    as no limit being taken.}
}
\value{
  a list containing some of the following components:
  \item{zeroth}{the value of the cumulant function at \eqn{\theta}{theta}.}
  \item{first}{the value of the first derivative at \eqn{\theta}{theta}, a
    vector having the same dimension as \eqn{\theta}{theta}.}
  \item{second}{the value of the second derivative at \eqn{\theta}{theta}, a
    \eqn{d \times d}{d by d} matrix if \eqn{d} is the dimension
    of \eqn{\theta}{theta} or a scalar if \eqn{\theta}{theta} is scalar.}
  \item{third}{the value of the third derivative at \eqn{\theta}{theta}, a
    \eqn{d \times d \times d}{d by d by d} array if \eqn{d} is the dimension
    of \eqn{\theta}{theta} or a scalar if \eqn{\theta}{theta} is scalar.}
}
\section{Note}{
Not intended for use by ordinary users.  Provides R interface for testing
to C code called by many other functions in the package.
}
\seealso{
\code{\link{families}}
}
\examples{
cumulant(-0.5, fam.bernoulli(), deriv = 3)
cumulant(-0.5, fam.bernoulli(), deriv = 3, delta = 1)
}
\keyword{misc}
