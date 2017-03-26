\name{families}
\alias{families}
\alias{fam.bernoulli}
\alias{fam.poisson}
\alias{fam.zero.truncated.poisson}
\alias{fam.normal.location.scale}
\alias{fam.multinomial}
\title{Families for Aster Models}
\usage{
fam.bernoulli()
fam.poisson()
fam.zero.truncated.poisson()
fam.normal.location.scale()
fam.multinomial(dimension)
}
\description{
Families known to the package.
These functions construct simple family specifications used
in specifying aster models.  Statistical properties of these families
are described.
}
\arguments{
  \item{dimension}{the dimension (number of categories) for the multinomial
    distribution.}
}
\details{
Currently implemented families are
  \describe{
    \item{\code{"bernoulli"}}{Bernoulli (binomial with sample size one).
      The distribution of any
      zero-or-one-valued random variable \eqn{Y}, which is the canonical
      statistic.  The mean value parameter is
      \deqn{\mu = E(Y) = \Pr(Y = 1).}{mu = E(Y) = Pr(Y = 1).}
      The canonical parameter is
      \eqn{\theta = \log(\mu) - \log(1 - \mu)}{theta = log(mu) - log(1 - mu)},
      also called logit of \eqn{\mu}{mu}.  The cumulant function is
      \deqn{c(\theta) = \log(1 + e^\theta).}{c(theta) = log(1 + exp(theta)).}

      This distribution has degenerate limiting distributions.  The lower
      limit as \eqn{\theta \to - \infty}{theta to minus infinity} is the
      distribution concentrated at zero, having cumulant function which
      is the constant function everywhere equal to zero.  The upper
      limit as \eqn{\theta \to + \infty}{theta to plus infinity} is the
      distribution concentrated at one, having cumulant function which
      is the identity function satisfying
      \eqn{c(\theta) = \theta}{c(theta) = theta} for all \eqn{\theta}{theta}.

      For predecessor (sample size) \eqn{n}, the successor is the sum of
      \eqn{n} independent and identically distributed (IID) Bernoulli
      random variables, that is,
      binomial with sample size \eqn{n}.  The mean value parameter is \eqn{n}
      times the mean value parameter for sample size one; the cumulant function
      is \eqn{n} times the cumulant function for sample size one; the canonical
      parameter is the same for all sample sizes.}
    \item{\code{"poisson"}}{Poisson.  The mean value parameter
      \eqn{\mu}{mu} is the mean of the Poisson distribution.
      The canonical parameter is \eqn{\theta = \log(\mu)}{theta = log(mu)}.
      The cumulant function is
      \deqn{c(\theta) = e^\theta.}{c(theta) = exp(theta).}

      This distribution has a degenerate limiting distribution.  The lower
      limit as \eqn{\theta \to - \infty}{theta to minus infinity} is the
      distribution concentrated at zero, having cumulant function which
      is the constant function everywhere equal to zero.  There is no upper
      limit because the canonical statistic is unbounded above.

      For predecessor (sample size) \eqn{n}, the successor is the sum of
      \eqn{n} IID Poisson random variables, that is,
      Poisson with mean \eqn{n \mu}{n mu}.  The mean value parameter is \eqn{n}
      times the mean value parameter for sample size one; the cumulant function
      is \eqn{n} times the cumulant function for sample size one; the canonical
      parameter is the same for all sample sizes.}
    \item{\code{"zero.truncated.poisson"}}{Poisson conditioned on being
      greater than zero.  Let \eqn{m} be the mean of the corresponding
      untruncated Poisson distribution.  Then the canonical parameters for both
      truncated and untruncated distributions are the same
      \eqn{\theta = \log(m)}{theta = log(m)}.
      The mean value parameter for the zero-truncated Poisson distribution is
      \deqn{\mu = \frac{m}{1 - e^{- m}}}{mu = 1 / (1 - exp(- m))}
      and the cumulant function is
      \deqn{c(\theta) = m + \log(1 - e^{- m}),}{m + log(1 - exp(- m)),}
      where \eqn{m} is as defined above,
      so \eqn{m = e^\theta}{m = exp(theta)}.

      This distribution has a degenerate limiting distribution.  The lower
      limit as \eqn{\theta \to - \infty}{theta to minus infinity} is the
      distribution concentrated at one, having cumulant function which
      is the identity function satisfying
      \eqn{c(\theta) = \theta}{c(theta) = theta} for all \eqn{\theta}{theta}.
      There is no upper
      limit because the canonical statistic is unbounded above.

      For predecessor (sample size) \eqn{n}, the successor is the sum of
      \eqn{n} IID zero-truncated Poisson random variables, which is not
      a brand-name distribution.  The mean value parameter is \eqn{n}
      times the mean value parameter for sample size one; the cumulant function
      is \eqn{n} times the cumulant function for sample size one; the canonical
      parameter is the same for all sample sizes.}
    \item{\code{"normal.location.scale"}}{The distribution of a normal
      random variable \eqn{X} with unknown mean \eqn{m} and unknown variance
      \eqn{v}.  Thought of as an exponential family, this is
      a two-parameter family, hence must have a two-dimensional canonical
      statistic \eqn{Y = (X, X^2)}.  The canonical parameter
      vector \eqn{\theta}{theta} has components
      \deqn{\theta_1 = \frac{m}{v}}{theta[1] = m / v}
      and
      \deqn{\theta_2 = - \frac{1}{2 v}.}{theta[2] = - 1 / (2 v).}
      The value of \eqn{\theta_1}{theta[1]} is unrestricted, but
      \eqn{\theta_2}{theta[2]} must be strictly negative.
      The mean value parameter vector \eqn{\mu}{mu} has components
      \deqn{\mu_1 = m = - \frac{\theta_1}{2 \theta_2}}{mu[1] = m = - theta[1] / (2 theta[2])}
      and
      \deqn{\mu_2 = v + m^2 = - \frac{1}{2 \theta_2} + \frac{\theta_1^2}{4 \theta_2^2}.}{mu[2] = v + m^2 = - 1 / (2 theta[2]) + theta[1]^2 / (4 theta[2]^2).}
      The cumulant function is
      \deqn{c(\theta) = - \frac{\theta_1^2}{4 \theta_2} + \frac{1}{2} \log\left(- \frac{1}{2 \theta_2}\right).}{c(theta) = - theta[1]^2 / (4 theta[2]) + (1 / 2) log(- 1 / (2 theta[2])).}

      This distribution has no degenerate limiting distributions, because
      the canonical statistic is a continuous random vector so the boundary
      of its support has probability zero.

      For predecessor (sample size) \eqn{n}, the successor is the sum of
      \eqn{n} IID random vectors \eqn{(X_i, X_i^2)}{(x[i], x[i]^2)},
      where each \eqn{X_i}{x[i]} is normal
      with mean \eqn{m} and variance \eqn{v}, and this is not
      a brand-name multivariate distribution (the first component of the sum
      is normal, the second component noncentral chi-square, and the
      components are not independent).
      The mean value parameter vector is \eqn{n}
      times the mean value parameter vector for sample size one;
      the cumulant function
      is \eqn{n} times the cumulant function for sample size one; the canonical
      parameter vector is the same for all sample sizes.}
    \item{\code{"multinomial"}}{Multinomial with sample size one.
      The distribution of any random vector \eqn{Y} having all components zero
      except for one component which is one (\eqn{Y} is the
      canonical statistic vector).
      The mean value parameter is the vector \eqn{\mu = E(Y)}{mu = E(Y)} having
      components
      \deqn{\mu_i = E(Y_i) = \Pr(Y_i = 1).}{mu[i] = E(Y[i]) = Pr(Y[i] = 1).}
      The mean value parameter vector \eqn{\mu}{mu} is given as a function
      of the canonical parameter vector \eqn{\theta}{theta} by
      \deqn{\mu_i = \frac{e^{\theta_i}}{\sum_{j = 1}^d e^{\theta_j}},}{mu[i] = exp(theta[i]) / sum(exp(theta)),}
      where \eqn{d} is the dimension of \eqn{Y} and \eqn{\theta}{theta}
      and \eqn{\mu}{mu}.  This transformation is not one-to-one;
      adding the same number
      to each component of \eqn{\theta}{theta} does not change the value
      of \eqn{\mu}{mu}.
      The cumulant function is
      \deqn{c(\theta) = \log\left(\sum_{j = 1}^d e^{\theta_j}\right).}{c(theta) = log(sum(exp(theta))).}

      This distribution is degenerate.  The sum of the components of the
      canonical statistic is equal to one with probability
      one, which implies the nonidentifiability of the \eqn{d}-dimensional
      canonical parameter vector mentioned above.  Hence one parameter
      (at least) is always constrained to to be zero in
      fitting an aster model with a multinomial family.

      This distribution has many degenerate distributions.  For any vector
      \eqn{\delta}{delta} the limit of distributions having canonical
      parameter vectors \eqn{\theta + s \delta}{theta + s delta} as
      \eqn{s \to \infty}{s to infinity} exists and is another
      multinomial distribution (the limit distribution in the direction
      \eqn{\delta}{delta}).
      Let \eqn{A} be the set of \eqn{i} such that 
      \eqn{\delta_i = \max(\delta)}{delta[i] = max(delta)},
      where \eqn{\max(\delta)}{max(delta)} denotes the maximum over the
      components of \eqn{\delta}{delta}.
      Then the limit distribution in the direction \eqn{\delta}{delta}
      has components \eqn{Y_i}{Y[i]} of the canonical statistic
      for \eqn{i \notin A}{i not in A} concentrated at zero.
      The cumulant function of this degenerate distribution is
      \deqn{c(\theta) = \log\left(\sum_{j \in A} e^{\theta_j}\right).}{c(theta) = log(sum(\{ theta[j] : j in A \})).}
      The canonical parameters \eqn{\theta_j}{theta[j]} for \eqn{j \notin A}{j not in A} are not identifiable, and one other canonical parameter is not
      identifiable because of the constraint that the sum of the components
      of the canonical statistic is equal to one with probability one.

      For predecessor (sample size) \eqn{n}, the successor is the sum of
      \eqn{n} IID multinomial-sample-size-one random vectors, that is,
      multinomial with sample size \eqn{n}.  The mean value parameter is \eqn{n}
      times the mean value parameter for sample size one; the cumulant function
      is \eqn{n} times the cumulant function for sample size one; the canonical
      parameter is the same for all sample sizes.}
  }
}
\value{
    a list of class \code{"astfam"} giving name and values of any
    hyperparameters.
}
\examples{
fam.bernoulli()
fam.multinomial(4)
}
\keyword{misc}
