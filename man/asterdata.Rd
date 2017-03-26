\name{asterdata}
\alias{asterdata}
\alias{validasterdata}
\alias{is.validasterdata}
\title{Object Describing Saturated Aster Model}
\description{
Functions to construct and test conformance to the contract for objects
of class \code{"asterdata"}.  All other functions in this package take
model descriptions of this form.
}
\usage{
asterdata(data, vars, pred, group, code, families, delta,
  response.name = "resp", varb.name = "varb",
  tolerance = 8 * .Machine$double.eps)
validasterdata(object, tolerance = 8 * .Machine$double.eps)
is.validasterdata(object, tolerance = 8 * .Machine$double.eps)
}
\arguments{
  \item{data}{a data frame containing response and predictor variables for
    the aster model.}
  \item{vars}{a character vector containing names of variables in the data
    frame \code{data} that are components of the response vector of the
    aster model.}
  \item{pred}{an integer vector satisfying \code{length(pred) == length(vars)}
    specifying the arrows of the subgraph of the aster model corresponding
    to a single individual.  Must be nonnegative and satisfy
    \code{all(pred < seq(along = pred))}.
    A zero value of \code{pred[j]} indicates the predecessor of node \code{j}
    is an initial node (formerly called root node) of the subgraph.
    A nonzero value of \code{pred[j]} indicates the predecessor of node
    \code{j} is node \code{pred[j]}.  In either case there is an arrow in
    the subgraph from predecessor node to successor node.}
  \item{group}{an integer vector satisfying \code{length(group) == length(vars)}
    specifying the lines of the subgraph of the aster model corresponding to
    a single individual, which in turn specify the dependence groups.
    Must be nonnegative and satisfy \code{all(group < seq(along = group))}.
    Nonzero elements of \code{group} indicate nodes of the subgraph that
    are connected by a line and hence are in the same dependence group:
    nodes \code{j} and \code{group[j]} are connected by a line.  Since
    nodes in the same dependence group must have the same predecessor,
    this requires \code{pred[group[j]] == pred[j]}.  Since
    nodes in the same dependence group must be in the same family,
    this requires \code{code[group[j]] == code[j]}.
    It also requires that the dimension of the family specified by
    \code{code[j]} be the same as the number of nodes in the dependence
    group.  Zero elements of \code{group} indicate nothing about dependence
    groups.

    The lines indicate a transitive relation.  If there is a line from
    node \code{j1} to node \code{j2} and a line from
    node \code{j2} to node \code{j3} then there is also a line from
    node \code{j1} to node \code{j3}, but this line need not be specified
    by the \code{group} vector, and indeed cannot.  If there is a dependence
    group with \code{d} nodes, then there are \code{choose(d, 2)} lines
    connecting these nodes, but the \code{group} vector can only specify
    \code{d - 1} lines which imply the rest.

    For example, if nodes \code{j1}, \code{j2}, \code{j3}, and \code{j4}
    are to make up a four-dimensional dependence group and \code{j1 < j2},
    \code{j2 < j3}, and \code{j3 < j4}, we must have \code{group[j1] == 0},
    \code{group[j2] == j1}, \code{group[j3] == j2}, and
    \code{group[j4] == j3}.
    This is forced by the requirement \code{all(group < seq(along = group))}.}
  \item{code}{an integer vector satisfying \code{length(code) == length(vars)}
    specifying the families corresponding to the dependence groups.
    This requires \preformatted{all(code \%in\% seq(along = families)}
    Node \code{j} is in a dependence group
    with family described by \code{families[code[j]]}.

    Note that \code{group[j] == k}
    requires \code{families[j] == families[k]} when \code{k != 0}.}
  \item{families}{a list of family specifications
    (see \code{\link{families}}).  Specifications of families not having
    hyperparameters may be abbreviated as character strings, for example,
    \code{"binomial"} rather than \code{fam.binomial()}.}
  \item{delta}{a numeric vector satisfying \code{length(delta) == length(vars)}
    specifying the degeneracies of the aster model for a single individual.
    The model specified is the limit as \eqn{s \to \infty}{s to infinity} of
    nondegenerate models having conditional canonical parameter vector
    \eqn{\theta + s \delta}{theta + s delta} (note that the conditional
    canonical parameter vector is always used here, regardless of whether
    conditional or unconditional canonical affine submodels are to be used).
    May be missing (and usually is) in which case \eqn{\delta = 0}{delta = 0}
    is implied, meaning the limit is trivial (same as not taking a limit).}
  \item{response.name}{a character string giving the name of the response
    vector.}
  \item{varb.name}{a character string giving the name of the factor covariate
    that says which of the variables in the data frame \code{data} correspond
    to which components of the response vector.}
  \item{tolerance}{numeric >= 0.  Relative errors smaller
    than \code{tolerance} are not considered in checking validity
    of normal location-scale data.}
  \item{object}{an object of class \code{"asterdata"}.  The function
    \code{validasterdata} always returns \code{TRUE} or throws an error with
    an informative message.  The function \code{is.validasterdata} never throws
    an error unless \code{object} has the wrong class, returning \code{TRUE}
    or \code{FALSE} according to whether \code{object} does or does not
    conform to the contract for class \code{"asterdata"}.}
}
\details{
  Response variables in dependence groups are taken to be in the order they
  appear in the response vector.  The first to appear in the response vector
  is the first canonical statistic for the dependence group distribution,
  the second to appear the second canonical statistic, and so forth.
  The number of response variables in the dependence group must match the
  dimension of the dependence group distribution.

  This function only handles the usual case where the subgraph for every
  individual is isomorphic to subgraph for every other individual
  and all initial nodes (formerly
  called root nodes) correspond to the constant one.  Each row of \code{data}
  is the data for one individual.  The vectors \code{vars}, \code{pred},
  \code{group}, \code{code}, and \code{delta} (if not missing) describe
  the subgraph for one individual (which is the same for all individuals).

  In other cases for which this function does not have the flexibility to
  construct the appropriate object of class \code{"asterdata"}, such an
  object will have to be constructed \dQuote{by hand} using R statements
  not involving this function or modifying an object produced by this
  function.  See the following section for description of such objects.
  The functions \code{validasterdata} and \code{is.validasterdata} can be
  used to check whether objects constructed \dQuote{by hand} have been
  constructed correctly.
}
\value{
  an object of class \code{"asterdata"} is a list containing the
  following components
    \item{redata}{a data frame having \code{nrow(data) * length(vars)} rows
      and containing variables having names
      in \code{setdiff(names(data), vars)} and also the names
      \code{"id"}, \code{response.name}, and \code{varb.name}.
      Produced from \code{data} using the \code{\link{reshape}}
      function.  Each variable in \code{setdiff(names(data), vars)} is repeated
      \code{length(vars)} times.  The variable named \code{response.name}
      is the concatenation of the variables in \code{data} with names
      in \code{vars}.  The variable named \code{varb.name} is a factor
      having levels \code{vars} that says which of the variables in the data
      frame \code{data} correspond to which components of the response vector.
      The variable named \code{"id"} is an integer vector that says which of
      the individuals (which rows of \code{data}) correspond to which rows
      of \code{redata}.  Not all objects of class \code{"asterdata"} need
      have an \code{id} variable, although all those constructed by this
      function do.}
    \item{repred}{an integer vector satisfying
      \code{length(repred) == nrow(redata)} specifying the arrows of the
      graph of the aster model for all individuals.
      Must be nonnegative and satisfy \code{all(repred < seq(along = repred))}.
      A zero value of \code{repred[j]} indicates the predecessor of node
      \code{j} is an initial node (formerly called root node) of the graph.
      A nonzero value of \code{repred[j]} indicates the predecessor of node
      \code{j} is node \code{repred[j]}.  In either case there is an arrow in
      the graph from predecessor node to successor node.

      Note that
      \code{repred} is determined by \code{pred} but is quite different from
      it.  Firstly, the lengths differ.  Secondly, \code{repred} is not just
      a repetition of \code{pred}.  The numbers in \code{pred}, if nonzero,
      are indices for the vector \code{vars} whereas the numbers
      in \code{repred}, if nonzero,
      are row indices for the data frame \code{redata}.}
    \item{initial}{a numeric vector specifying constants associated with
      initial nodes (formerly called root nodes) of the graphical model
      for all individuals.  If \code{repred[j] == 0} then the predecessor
      of node \code{j} is an initial node associated with the constant
      \code{initial[j]}, which must be a positive integer unless the
      family associated with the arrow from this initial node to node \code{j}
      is infinitely divisible (the only such family currently implemented
      being Poisson), in which case \code{initial[j]} must
      be a strictly positive and finite real number.  If \code{repred[j] != 0},
      then \code{initial[j]} is ignored and may be any numeric value, including
      \code{NA} or \code{NaN}.  This function always makes \code{initial}
      equal to \code{rep(1, nrow(redata))} but the more general description
      above is valid for objects of class \code{"asterdata"} constructed
      \dQuote{by hand}.}
    \item{regroup}{an integer vector satisfying
      \code{length(regroup) == nrow(redata)}
      specifying the lines of the graph of the aster model for all individuals,
      which in turn specify the dependence groups.
      Must be nonnegative
      and satisfy \code{all(regroup < seq(along = regroup))}.
      Nonzero elements of \code{regroup} indicate nodes of the graph that
      are connected by a line and hence are in the same dependence group:
      nodes \code{j} and \code{regroup[j]} are connected by a line.  Since
      nodes in the same dependence group must have the same predecessor,
      this requires \code{repred[regroup[j]] == repred[j]}.  Since
      nodes in the same dependence group must be in the same family,
      this requires \code{recode[regroup[j]] == recode[j]}.

      It also requires that the dimension of the family specified by
      \code{recode[j]} be the same as the number of nodes in the dependence
      group.  Zero elements of \code{regroup} indicate nothing about dependence
      groups.

      The lines indicate a transitive relation.  If there is a line from
      node \code{j1} to node \code{j2} and a line from
      node \code{j2} to node \code{j3} then there is also a line from
      node \code{j1} to node \code{j3}, but this line need not be specified
      by the \code{group} vector, and indeed cannot.  If there is a dependence
      group with \code{d} nodes, then there are \code{choose(d, 2)} lines
      connecting these nodes, but the \code{group} vector can only specify
      \code{d - 1} lines which imply the rest.

      For example, if nodes \code{j1}, \code{j2}, \code{j3}, and \code{j4}
      are to make up a four-dimensional dependence group and \code{j1 < j2},
      \code{j2 < j3}, and \code{j3 < j4}, we must have \code{regroup[j1] == 0},
      \code{regroup[j2] == j1}, \code{regroup[j3] == j2}, and
      \code{regroup[j4] == j3}.
      This is forced by the requirement
      \code{all(regroup < seq(along = regroup))}.

      Note that
      \code{regroup} is determined by \code{group} but is quite different from
      it.  Firstly, the lengths differ.  Secondly, \code{regroup} is not just
      a repetition of \code{group}.  The numbers in \code{group}, if nonzero,
      are indices for the vector \code{vars} whereas the numbers
      in \code{regroup}, if nonzero,
      are row indices for the data frame \code{redata}.}
    \item{recode}{an integer vector satisfying
      \code{length(recode) == nrow(redata)}
      specifying the families corresponding to the dependence groups.
      This requires \preformatted{all(recode \%in\% seq(along = families)}
      Node \code{j} is in a dependence group
      with family described by \code{families[recode[j]]}.

      Note that \code{regroup[j] == k}
      requires \code{recode[j] == recode[k]} when \code{regroup[j] != 0}.
      Also note that
      \code{recode} is determined by \code{code} but is different from
      it.  Firstly, the lengths differ.  Secondly, \code{recode} need not be
      just a repetition of \code{code}.
      This function always makes \code{recode}
      equal to \code{rep(code, each = nrow(redata))} but the more general
      description
      above is valid for objects of class \code{"asterdata"} constructed
      \dQuote{by hand}.}
    \item{families}{a copy of the argument of the same name of this function
      except that any character string abbreviations are converted to objects
      of class \code{"astfam"}.}
    \item{redelta}{a numeric vector satisfying
      \code{length(redelta) == nrow(redata)}
      specifying the degeneracies of the aster model for all individuals.
      If not the zero vector, the degenerate model
      specified is the limit as \eqn{s \to \infty}{s to infinity} of
      nondegenerate models having conditional canonical parameter vector
      \eqn{\theta + s \delta}{theta + s delta} (note that the conditional
      canonical parameter vector is always used here, regardless of whether
      conditional or unconditional canonical affine submodels are to be used).

      Note that
      \code{redelta} is determined by \code{delta} but is different from
      it.  Firstly, the lengths differ.  Secondly, \code{redelta} need not be
      just a repetition of \code{delta}.
      This function always makes \code{redelta}
      equal to \code{rep(delta, each = nrow(redata))} but the more general
      description
      above is valid for objects of class \code{"asterdata"} constructed
      \dQuote{by hand}.}
    \item{response.name}{a character string giving the name of the response
      variable in \code{redata}.  For this function, a copy of the argument
      \code{response.name}.}
    \item{varb.name}{a character string giving the name of the \dQuote{varb}
      variable in \code{redata}.  For this function, a copy of the argument
      \code{varb.name}.}
  In addition an object of class \code{"asterdata"} may contain (and those
  constructed by this function do contain) components
  \code{pred}, \code{group}, and \code{code},
  which are copies of the arguments of the same names of this function.
  Objects of class \code{"asterdata"} not constructed by this function need
  not contain these additional components, since they may make no sense if
  the graph for all individuals is not the repetition of isomorphic subgraphs,
  one for each individual.
}
\seealso{\code{\link{families}} and \code{\link{subset.asterdata}}}
\examples{
data(test1)
fred <- asterdata(test1, vars = c("m1", "n1", "n2"), pred = c(0, 1, 1),
    group = c(0, 0, 2), code = c(1, 2, 2),
    families = list("bernoulli", "normal.location.scale"))
is.validasterdata(fred)
}
\keyword{manip}
