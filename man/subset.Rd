\name{subset.asterdata}
\alias{subset.asterdata}
\title{Subset Object Describing Saturated Aster Model}
\description{
Subset an object of class \code{"asterdata"},
for which see \code{\link{asterdata}}.
}
\usage{
\method{subset}{asterdata}(x, subset, successors = TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"asterdata"}, for which see
    \code{\link{asterdata}}.}
  \item{subset}{a logical vector indicating nodes of the graph to keep:
    missing values are taken as false.}
  \item{successors}{a logical scalar indicating whether the subgraph must be
    a union of connected components of the original graph, that is, if all
    successors of nodes in the subset must also be in the subset.}
  \item{\dots}{further arguments, which are ignored (this argument is required
    for methods of the generic function \code{\link{subset}} but is not used
    for this method.)}
}
\details{
  Argument \code{subset} is a logical vector of the same length as the number
  of nodes in the graph specified by argument \code{x}.  It indicates the
  subset of nodes in the subgraph wanted.  The subgraph must be closed with
  respect to predecessors (all predecessors of nodes in the subset are also
  in the subset) and if \code{successors = TRUE} with respect
  to successors (all successors of nodes in the subset are also in the subset).
  And similarly for dependence groups: each dependence
  group in the original graph must have all or none of its elements
  in the subgraph.
}
\value{
  an object of class \code{"asterdata"} that represents the aster model having
  subgraph with nodes specified by \code{subset}.
}
\seealso{\code{\link{asterdata}}}
\examples{
data(echinacea)
#### select one individual from each level of pop
foo <- echinacea$redata$pop
bar <- match(levels(foo), as.character(foo))
baz <- is.element(echinacea$redata$id, echinacea$redata$id[bar])
out <- subset(echinacea, baz)
}
\keyword{manip}
