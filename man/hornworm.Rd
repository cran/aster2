
\name{hornworm}
\docType{data}
\alias{hornworm}
\title{Life History Data on Manduca sexta}
\description{
  Data on life history traits for the tobacco hornworm
  \emph{Manduca sexta}
}
\usage{data(hornworm)}
\format{
  An object of class \code{"asterdata"} (see \code{\link{asterdata}})
  comprising records for 162 insects (54 female, 68 male, and 40 for which
  there was no opportunity to determine sex) observed over 40 days.
  Nodes of the graph for one individual are associated with the variables
  (levels of the factor \code{hornworm$redata$varb}) in dependence groups
  \describe{
    \item{P}{Bernoulli.  Predecessor 1 (initial node).  Indicator of pupation.}
    \item{T330, T331, T332}{Three-dimensional multinomial dependence group.
      Predecessor \code{P}.}
      \describe{
      \item{T330}{Indicator of death after pupation.  In these data,
        all deaths after pupation are considered to have happened on
        day 33 regardless of when they occurred (because the actual
        day of death was not recorded in the original data).}
      \item{T331}{Indicator of survival to day 33 but still pre-eclosion.}
      \item{T332}{Indicator of eclosion (emergence from pupa as adult moth
         on day 33.}
      }
    \item{B33}{Zero-truncated Poisson.  Predecessor \code{T332}.
       Count of ovarioles on day 33.  Only females have this node in their
       graphs.}
    \item{Tx1, Tx2}{For \code{x} = 34, \ldots, 40. Two-dimensional multinomial
      dependence group.  Predecessor \code{Tw1}, where \code{w = x - 1}.}
      \describe{
        \item{Tx1}{Indicator of survival to day \code{x} but still
          pre-eclosion.}
        \item{Tx2}{Indicator of eclosion (emergence from pupa as adult moth
          on day \code{x}.}
      }
    \item{Bx}{Zero-truncated Poisson.  Predecessor \code{Tx2}.
       Count of ovarioles on day \code{x}.  Only females have these nodes in
       their graph.}
   }
   Covariates are
   \describe{
    \item{Sex}{a factor.  \code{F} is known female, \code{M} is known male,
      \code{U} is unknown (no opportunity to observe).}
    \item{Time_2nd}{time (in weeks) to reach the 2nd instar stage.  Larval
      instars are stages between molts (shedding of exoskeleton) of the larval
      form (caterpillar).}
    \item{Mass_2nd}{mass (in grams) at the 2nd instar stage.}
    \item{Mass_Repro}{mass (in grams) at eclosion.}
    \item{LarvaID}{name of an individual in the original data.}
  }
}
\details{
This is the data described by and analyzed by non-aster methods by
Kingsolver et al. (2012) and re-analyzed using this package by Eck et al.
(submitted).

For an illustration of the graph, see Figure 1 in Eck et al. (submitted).

In the description above, a concrete example of the \code{x} and \code{w}
notation is that T351 and T352 form a two-dimensional multinomial dependence
group, the predecessor of which is T341, and B35 is a dependence group all
by itself, its predecessor being T352.

Every multinomial dependence group acts like a switch.  If the predecessor
is one, the dependence group is multinomial with sample size one (exactly
one variable is one and the rest are zero).  So this indicates which way
the life history goes.  If the predecessor is zero, then all successors are
zero.  This goes for all variables in any aster model.  If \code{Tx2} is zero,
then so is \code{Bx}.  The ovariole count is zero except for the day on
which the individual eclosed.
}
\source{
Joel Kingsolver
\url{https://bio.unc.edu/people/faculty/kingsolver/}
}
\references{
Kingsolver, J. G., Diamond, S. E., Seiter, S. A., and Higgins, J. K. (2012)
Direct and indirect phenotypic selection on developmental trajectories
  in \emph{Manduca sexta}.
\emph{Functional Ecology} \bold{26} 598--607.

Eck, D., Shaw, R. G., Geyer, C. J., and Kingsolver, J. (submitted)
An integrated analysis of phenotypic selection on insect body size
   and development time.
}
\examples{
data(hornworm)
names(hornworm)
names(hornworm$redata)
levels(hornworm$redata$varb)
}
\keyword{datasets}

