\name{aster2-package}
\alias{aster2-package}
\alias{aster2}
\docType{package}
\title{
Aster Models
}
\description{
Aster models are exponential family graphical models that
combine aspects of generalized linear models and survival analysis.
}
\details{
\tabular{ll}{
Package: \tab aster2\cr
Version: \tab 0.1\cr
Date: \tab 2010-09-12\cr
Depends: \tab R (>= 2.6.0)\cr
Suggests: \tab aster\cr
License: \tab GPL (>= 2)\cr
URL: \tab \url{http://www.stat.umn.edu/geyer/aster/}\cr
}
This package is still under development, only about half finished.
However, it does do maximum likelihood for unconditional aster models
with dependence groups, which the old package \code{aster} does not.

The main differences between this package and the old package are as follows.
\enumerate{
\item The old package had triple indices for model matrices.  The first index
ran over individuals, the second index over nodes of the graph for
an individual, and the third index over regression coefficients.  Consequently
the model matrix was represented (sometimes, but not consistently) as a
three-dimensional array rather than a matrix, which was very confusing, even
to the package author.  This package ignores individuals, one index runs
over all nodes of the combined graph for all individuals.  Thus model matrices
are always matrices.
\item The old package did not implement dependence groups, although they were
described in Geyer, Wagenius and Shaw (2007).  This package does.
Consequently, this package requires a data frame, a vector \code{pred} that
indicates predecessors, a vector \code{group} that indicates individuals in
the same dependence group, and a vector \code{fam} that indicates families
to specify a saturated aster model (the old package required only
the data frame, \code{pred}, and \code{fam}).
To facilitate the old style model specification, there is a new function
\code{\link{asterdata}} that constructs objects of class \code{"asterdata"}
given an old style data frame, \code{pred}, and \code{fam}.  All other
functions of the package take objects of class \code{"asterdata"} as model
specifications.
\item The function \code{\link[aster]{predict.aster}} in the old package
did some parameter transformations, but not all, and the returned value,
when a list, had a component \code{gradient}, that was undocumented but useful
in applying the delta method.  The functions
\code{\link{transformSaturated}},
\code{\link{transformConditional}}, and
\code{\link{transformUnconditional}} in this package transform
between any of the following parameter vectors:
the conditional canonical parameter \eqn{\theta}{theta},
the unconditional canonical parameter \eqn{\varphi}{phi},
the conditional mean value parameter \eqn{\xi}{xi},
the unconditional mean value parameter \eqn{\mu}{mu},
the canonical affine submodel canonical parameter \eqn{\beta}{beta},
and (unconditional aster models only)
the canonical affine submodel mean value parameter \eqn{\tau}{tau}
(this last parameter is new, not discussed in the cited papers below, it is
\eqn{\tau = M^T \mu}{tau = M^T mu}, where \eqn{M} is the model matrix).
The change of parameter from \eqn{\tau}{tau} to \eqn{\beta}{beta}
is equivalent to maximum likelihood estimation for an unconditional
aster model when the value \eqn{\tau = M^T y}{tau = M^T y} is used,
where \eqn{y} is the response vector.  All of these transformation functions
also compute derivatives, if requested.  See examples.
}
}
\section{Bugs}{
Functions analogous to \code{aster}, \code{anova}, and \code{predict}
in the old package are missing, thus model fitting, hypothesis tests,
and confidence intervals are more cumbersome.  In fact, since there is
no function to calculate log likelihoods (like \code{mlogl} in the old
package), there is no way to do likelihood ratio tests (but Rao or Wald
tests could be done, for unconditional aster models, since the derivative
of the log likelihood is observed minus expected
\eqn{M^T (y - \mu)}{M^T (y - mu)}.
}
\references{
Geyer, C. J., Wagenius, S., and Shaw, R. G. (2007)
Aster Models for Life History Analysis.
\emph{Biometrika} \bold{94} 415--426.

Shaw, R. G., Geyer, C. J., Wagenius, S., Hangelbroek, H. H.
    and Etterson, J. R. (2008)
Unifying Life History Analyses for Inference of Fitness
    and Population Growth.
\emph{American Naturalist}, \bold{172}, E35--E47. 
}
\keyword{package}
\seealso{
\code{\link{asterdata}}, \code{\link{transformSaturated}},
\code{\link{families}}
}
\examples{
data(echinacea)
#### estimate MLE (simpler model than in Biometrika paper cited, not as good)
hdct <- as.numeric(grepl("hdct", as.character(echinacea$redata$varb)))
modmat <- model.matrix(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
    data = echinacea$redata)
tau.hat <- as.numeric(t(modmat) \%*\% echinacea$redata$resp)
beta.hat <- transformUnconditional(tau.hat, modmat, echinacea,
    from = "tau", to = "beta")
inverse.fisher <- jacobian(tau.hat, echinacea, transform = "unconditional",
    from = "tau", to = "beta", modmat = modmat)
#### now have MLE (beta.hat) and pseudo-inverse of Fisher information
#### (inverse.fisher), pseudo-inverse because modmat is not full rank
foo <- cbind(beta.hat, sqrt(diag(inverse.fisher)))
foo <- cbind(foo, foo[ , 1]/foo[ , 2])
foo <- cbind(foo, 2 * pnorm(- abs(foo[ , 3])))
dimnames(foo) <- list(colnames(modmat),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
printCoefmat(foo)
#### coeffients constrained to be zero because parameterization is not
#### identifiable have estimate zero and std. error zero (and rest NA)

#### estimate fitness in populations
#### generate new data with one individual in each pop at location (0, 0)
pop.names <- levels(echinacea$redata$pop)
pop.idx <- match(pop.names, as.character(echinacea$redata$pop))
pop.id <- echinacea$redata$id[pop.idx]
newdata <- subset(echinacea, echinacea$redata$id \%in\% pop.id)
newdata$redata[ , "nsloc"] <- 0
newdata$redata[ , "ewloc"] <- 0
hdct <- as.integer(grepl("hdct", as.character(newdata$redata$varb)))
#### modmat for new data
newmodmat <- model.matrix(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
    data = newdata$redata)
#### matrix that when multiplied mean value parameter vector gives fitness
#### in each pop
amat <- matrix(NA, nrow = length(pop.id), ncol = nrow(newmodmat))
for (i in 1:nrow(amat))
    amat[i, ] <- as.numeric(grepl(paste("^", pop.id[i], ".hdct", sep = ""),
        rownames(newmodmat)))
#### transform to expected fitness parameters
efit <- transformUnconditional(beta.hat, newmodmat, newdata,
    from = "beta", to = "mu")
efit <- as.numeric(amat \%*\% efit)
#### jacobian matrix of this transformation
jack <- jacobian(beta.hat, newdata, transform = "unconditional",
    from = "beta", to = "mu", modmat = newmodmat)
#### delta method standard errors
sefit <- sqrt(diag(amat \%*\% jack \%*\% inverse.fisher \%*\% t(jack) \%*\% t(amat)))
foo <- cbind(efit, sefit)
dimnames(foo) <- list(pop.names, c("Est. fitness", "Std. Error"))
print(foo)
}
