\name{test1}
\alias{test1}
\docType{data}
\title{Test Data}
\description{
Test data of no biological interest.  Does have all families
implemented at the time the test data was created.  No predictor variables.
}
\usage{data(test1)}
\format{
  A data frame with 100 observations on the following 8 variables.
  \describe{
    \item{\code{m1}}{a numeric vector, part of a multinomial dependence group
      (with \code{m2} and \code{m3}).  Predecessor of this group is the
      constant 1.}
    \item{\code{m2}}{a numeric vector.}
    \item{\code{m3}}{a numeric vector.}
    \item{\code{n1}}{a numeric vector, part of a normal location-scale
      dependence group (with \code{n2}).  Predecessor of this group
      is \code{m1}.}
    \item{\code{n2}}{a numeric vector (actually \code{n1^2}).}
    \item{\code{b1}}{a numeric vector, Bernoulli.  Predecessor is \code{m2}.}
    \item{\code{p1}}{a numeric vector, Poisson.  Predecessor is \code{m3}.}
    \item{\code{z1}}{a numeric vector, zero-truncated Poisson.  Predecessor
      is \code{b1}.}
  }
}
\source{
created by R script \code{test1.R} in directory \code{makedata} of the
installation directory for this package.
}
\examples{
data(test1)
fred <- asterdata(test1,
    vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
    pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
    code = c(1, 1, 1, 2, 2, 3, 4, 5),
    families = list(fam.multinomial(3), "normal.location.scale",
    "bernoulli", "poisson", "zero.truncated.poisson"))
}
\keyword{datasets}
