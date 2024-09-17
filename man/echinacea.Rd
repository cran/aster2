
\name{echinacea}
\docType{data}
\alias{echinacea}
\title{Life History Data on Echinacea angustifolia}
\description{
  Data on life history traits for the purple coneflower
  \emph{Echinacea angustifolia}
}
\usage{data(echinacea)}
\format{
  An object of class \code{"asterdata"} (see \code{\link{asterdata}})
  comprising  records for 570 plants observed over three years.
  Nodes of the graph for one individual are associated with the variables
  (levels of the factor \code{echinacea$redata$varb})
   \describe{
    \item{ld02}{Indicator of being alive in 2002.  Bernoulli,
      predecessor the constant one.}
    \item{ld03}{Ditto for 2003.  Bernoulli, predecessor \code{ld02}.}
    \item{ld04}{Ditto for 2004.  Bernoulli, predecessor \code{ld03}.}
    \item{fl02}{Indicator of flowering 2002.  Bernoulli,
      predecessor \code{ld02}.}
    \item{fl03}{Ditto for 2003.  Bernoulli, predecessor \code{ld03}.}
    \item{fl04}{Ditto for 2004.  Bernoulli, predecessor \code{ld04}.}
    \item{hdct02}{Count of number of flower heads in 2002.
      Zero-truncated Poisson, predecessor \code{fl02}.}
    \item{hdct03}{Ditto for 2003.
      Zero-truncated Poisson, predecessor \code{fl03}.}
    \item{hdct04}{Ditto for 2004.
      Zero-truncated Poisson, predecessor \code{fl04}.}
   }
   Covariates are
   \describe{
    \item{pop}{the remnant population of origin of the plant
       (all plants were grown together, \code{pop} encodes ancestry).}
    \item{ewloc}{east-west location in plot.}
    \item{nsloc}{north-south location in plot.}
  }
}
\details{
This is the data for the example in Geyer, Wagenius, and Shaw (2007).
These data were included in the R package \code{aster} which was the
predecessor of this package as the dataset \code{echinacea}.
}
\source{
Stuart Wagenius,
\url{https://www.chicagobotanic.org/research/staff/wagenius}
}
\references{
Geyer, C. J., Wagenius, S., and Shaw, R. G. (2007)
Aster Models for Life History Analysis.
\emph{Biometrika} \bold{94} 415--426.
}
\examples{
data(echinacea)
names(echinacea)
names(echinacea$redata)
levels(echinacea$redata$varb)
}
\keyword{datasets}

