\name{Transform}
\alias{Transform}
\alias{transformSaturated}
\alias{transformConditional}
\alias{transformUnconditional}
\alias{jacobian}
\alias{validtheta}
\alias{is.validtheta}
\alias{validxi}
\alias{is.validxi}
\title{Change-of-Parameter Functions for Aster Models}
\usage{
transformSaturated(parm, data, from = c("theta", "phi", "xi", "mu"),
  to = c("theta", "phi", "xi", "mu"), differential,
  tolerance = 8 * .Machine$double.eps)
transformConditional(parm, modmat, data, from = "beta",
  to = c("theta", "phi", "xi", "mu"), differential,
  offset, tolerance = 8 * .Machine$double.eps)
transformUnconditional(parm, modmat, data, from = c("beta", "tau"),
  to = c("beta", "theta", "phi", "xi", "mu", "tau"),
  differential, offset, tolerance = 8 * .Machine$double.eps)
jacobian(parm, data,
    transform = c("saturated", "conditional", "unconditional"),
    from = c("beta", "theta", "phi", "xi", "mu", "tau"),
    to = c("beta", "theta", "phi", "xi", "mu", "tau"),
    modmat, offset, tolerance = 8 * .Machine$double.eps)
validtheta(data, theta, tolerance = 8 * .Machine$double.eps)
is.validtheta(data, theta, tolerance = 8 * .Machine$double.eps)
validxi(data, xi, tolerance = 8 * .Machine$double.eps)
is.validxi(data, xi, tolerance = 8 * .Machine$double.eps)
}
\description{
Calculate a change-of-parameter for a saturated aster model or the
derivative of such a change-of-parameter.  Validate certain parameter vectors.
}
\arguments{
  \item{parm}{parameter vector to transform,
    a numerical vector of length \code{nrow(data$redata)}
    for \code{transformSaturated} or of length \code{ncol(modmat)}
    for \code{transformConditional} and \code{transformUnconditional}.}
  \item{data}{an object of class \code{"asterdata"} produced by
    \code{\link{asterdata}} or \dQuote{by hand} such that
    \code{is.validasterdata(data)} returns \code{TRUE}.  The specification
    of the aster model.}
  \item{from}{the kind of parameter which \code{parm} is.  May be abbreviated.}
  \item{to}{the kind of parameter to which \code{parm} is to be converted.
    May be abbreviated.}
  \item{differential}{if not missing a numeric vector of the same length
    as \code{parm}.  If missing the change-of-parameter is calculated.
    If not missing the derivative of the change-of-parameter is calculated.}
  \item{modmat}{the model matrix for a canonical affine submodel, a
    numerical matrix having \code{nrow(data$redata)} rows
    and \code{length(beta)} columns for \code{transformConditional}
    or \code{length(parm)} columns for \code{transformUnconditional}.}
  \item{offset}{the offset vector for a canonical affine submodel, a
    a numerical vector of length \code{nrow(data$redata)}.  May be missing,
    in which case offset vector equal to zero is used.}
  \item{theta}{conditional canonical parameter vector to validate,
    a numerical vector of length \code{nrow(data$redata)}.}
  \item{xi}{conditional canonical parameter vector to validate,
    a numerical vector of length \code{nrow(data$redata)}.}
  \item{tolerance}{numeric >= 0.  Relative errors smaller
    than \code{tolerance} are not considered in checking validity
    of \code{xi} for multinomial data.}
  \item{transform}{the \dQuote{transform} function that will be called to
    calculate derivatives, e. g., \code{transform == "saturated"} means the
    function \code{transformSaturated} will be called.  May be abbreviated.}
}
\details{
If \code{differential} is missing, the returned value is a new parameter
vector of the specified type.  If \code{differential} is not missing,
the returned value is the derivative evaluated at \code{parm}
and \code{differential}, that is, if \eqn{f} is the change-of variable
function and \eqn{\psi}{psi} is the \code{from} parameter, then
\eqn{f(\psi)}{f(psi)} is calculated when the differential is missing and
\eqn{f'(\psi)(\delta)}{f'(psi)(delta)} is calculated when the
differential \eqn{\delta}{delta} is not missing, where the latter is defined by
\deqn{f(\psi + \delta) \approx f(\psi) + f'(\psi)(\delta)}{f(psi + delta) approx f(psi) + f'(psi)(delta)} for small \eqn{\delta}{delta}.

The kinds of parameters are \code{"theta"} the conditional canonical parameter
for the saturated model, \code{"phi"} the unconditional canonical parameter
for the saturated model, \code{"xi"} the conditional mean value parameter
for the saturated model, \code{"mu"} the unconditional mean value parameter
for the saturated model,
\code{"beta"} the regression coefficient parameter for a canonical affine
submodel (\eqn{\theta = a + M \beta}{theta = a + M beta} for a conditional
canonical affine submodel or
\eqn{\varphi = a + M \beta}{phi = a + M beta} for an unconditional
canonical affine submodel, where \eqn{a} is the offset vector
and \eqn{M} is the model matrix),
\code{"tau"} the mean value parameter for an unconditional canonical affine
submodel (\eqn{\tau = M^T \mu}{tau = M^T mu},
where \eqn{M} is the model matrix).

Only the conditional canonical parameter vector \eqn{\theta}{theta} and
the conditional mean value parameter vector \eqn{\xi}{xi} can be checked
directly.  (To check the validity of another parameter, transform to one
of these and check that.)  This means that in conversions to these parameters
the output vector is checked rather than the input vector, and conversions
(apparently) not involving these parameters (which do go through these
parameters inside the transformation function) a conversion to one of
these parameters is what is checked rather than the input vector.

The function \code{jacobian} turns the derivative considered as
a linear transformation calculated by the \dQuote{transform} functions
into the matrix that represents the linear transformation (sometimes
called the Jacobian matrix of the transformation).  The arguments
\code{modmat} and \code{offset} are only used if
\code{transform == "conditional"} or \code{transform == "unconditional"},
and (as with the \dQuote{transform} functions) the argument \code{offset}
may be missing, in which case the zero vector is used.  Not all of the
candidate values for \code{from} and \code{to} arguments
for the \code{jacobian} function are valid: the value must be valid for
the \dQuote{transform} function that will be called.
}
\value{
  a numeric vector of the same length as \code{parm}.  The new parameter if
  \code{deriv == FALSE} or the transform of the differential
  if \code{deriv = TRUE}.  See details.
}
\seealso{
\code{\link{asterdata}}
}
\examples{
data(echinacea)
theta <- rnorm(nrow(echinacea$redata), 0, 0.1)
phi <- transformSaturated(theta, echinacea, from = "theta", to = "phi")
## rarely (if ever) want jacobian for unsaturated model transform
## result here is 5130 by 5130 matrix
\dontrun{jack <- jacobian(theta, echinacea, from = "theta", to = "phi")}
}
\keyword{misc}
