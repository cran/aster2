\name{constancy}
\alias{constancy}
\alias{is.same}
\title{Constancy Spaces for Aster Models}
\usage{
constancy(data, parm.type = c("theta", "phi"))
is.same(parm1, parm2, data, parm.type = c("theta", "phi"),
    tolerance = sqrt(.Machine$double.eps))
}
\description{
Produce basis for constancy space of an aster model.
Test whether the difference of two canonical parameter vectors is
in the constancy space (so the two parameter vectors correspond
to the same probability model).
}
\arguments{
  \item{data}{an object of class \code{"asterdata"} produced by
    \code{\link{asterdata}} or \dQuote{by hand} such that
    \code{is.validasterdata(data)} returns \code{TRUE}.  The specification
    of the aster model.}
  \item{parm.type}{the parametrization for which the constancy space is wanted.}
  \item{parm1}{a parameter vector of the type specified by \code{parm.type}
    for the saturated aster model specified by \code{data}.}
  \item{parm2}{another parameter vector of the type specified
    by \code{parm.type} for the saturated aster model specified
    by \code{data}.}
  \item{tolerance}{numeric >= 0.  Relative errors smaller
    than \code{tolerance} are not considered in the comparison.}
}
\details{
There is no need for functions to test whether different mean value parameters
(\eqn{\xi}{xi} or \eqn{\mu}{mu}) correspond to the same probability
distribution because these parametrizations are identifiable (different valid
parameter vectors correspond to different probability distributions).
}
\value{
  for \code{is.same} a logical value;
  for \code{constancy}
  a matrix whose rows constitute a basis for the constancy space.
  This means that if \eqn{\delta}{delta} is a linear combination of rows
  of this matrix then for all real \eqn{s} the distributions having parameter
  vectors \eqn{\psi}{psi} and \eqn{\psi + s \delta}{psi + s * delta} are the
  same, where \eqn{\psi = \theta}{psi = theta}
  or \eqn{\psi = \varphi}{psi = phi} depending on whether
  \code{parm.type = "theta"} or \code{parm.type = "phi"}.
  Conversely, if \eqn{\psi_1}{psi1} and \eqn{\psi_2}{psi2} are valid parameter
  vectors of the same type, then they correspond to the same probability
  distribution only if \eqn{\psi_1 - \psi_2}{psi1 - psi2} is a linear
  combination of rows of this matrix.
}
\seealso{
\code{\link{asterdata}}
}
\examples{
data(test1)
fred <- asterdata(test1,
    vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
    pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
    code = c(1, 1, 1, 2, 2, 3, 4, 5),
    families = list(fam.multinomial(3), "normal.location.scale",
    "bernoulli", "poisson", "zero.truncated.poisson"))
cmat <- constancy(fred, parm.type = "phi")
}
\keyword{misc}
