
 library(aster2)

 ##### Bernoulli, lower limit #####

 theta <- seq(-3, 3, 0.1)

 cumfun <- function(theta) rep(0, length(theta))
 moofun <- cumfun
 voofun <- cumfun
 thoofun <- cumfun

 zeroth <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = -1)$zeroth, theta)
 all.equal(zeroth, cumfun(theta))

 first <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = -1, deriv = 1)$first, theta)
 all.equal(first, moofun(theta))

 second <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = - 1, deriv = 2)$second, theta)
 all.equal(second, voofun(theta))

 third <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = - 1, deriv = 3)$third, theta)
 all.equal(third, thoofun(theta))

 foo <- link(0, fam.bernoulli(), delta = -1, deriv = 1)
 # does not matter what zeroth is, theta is not identifiable for
 # degenerate model
 is.finite(foo$zeroth)
 all.equal(0, foo$first)

 ##### Bernoulli, upper limit #####

 cumfun <- function(theta) theta
 moofun <- function(theta) rep(1, length(theta))
 voofun <- function(theta) rep(0, length(theta))
 thoofun <- voofun

 zeroth <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = 1)$zeroth, theta)
 all.equal(zeroth, cumfun(theta))

 first <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = 1, deriv = 1)$first, theta)
 all.equal(first, moofun(theta))

 second <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = 1, deriv = 2)$second, theta)
 all.equal(second, voofun(theta))

 third <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     delta = 1, deriv = 3)$third, theta)
 all.equal(third, thoofun(theta))

 foo <- link(1, fam.bernoulli(), delta = 1, deriv = 1)
 # does not matter what zeroth is, theta is not identifiable for
 # degenerate model
 is.finite(foo$zeroth)
 all.equal(0, foo$first)

 ##### Poisson, lower limit #####

 theta <- seq(-3, 3, 0.1)

 cumfun <- function(theta) rep(0, length(theta))
 moofun <- cumfun
 voofun <- cumfun
 thoofun <- cumfun

 zeroth <- mapply(function(theta) cumulant(theta, fam.poisson(),
     delta = -1)$zeroth, theta)
 all.equal(zeroth, cumfun(theta))

 first <- mapply(function(theta) cumulant(theta, fam.poisson(),
     delta = -1, deriv = 1)$first, theta)
 all.equal(first, moofun(theta))

 second <- mapply(function(theta) cumulant(theta, fam.poisson(),
     delta = - 1, deriv = 2)$second, theta)
 all.equal(second, voofun(theta))

 third <- mapply(function(theta) cumulant(theta, fam.poisson(),
     delta = - 1, deriv = 3)$third, theta)
 all.equal(third, thoofun(theta))

 foo <- link(0, fam.poisson(), delta = -1, deriv = 1)
 # does not matter what zeroth is, theta is not identifiable for
 # degenerate model
 is.finite(foo$zeroth)
 all.equal(0, foo$first)

 ##### zero-truncated Poisson, lower limit #####

 cumfun <- function(theta) theta
 moofun <- function(theta) rep(1, length(theta))
 voofun <- function(theta) rep(0, length(theta))
 thoofun <- voofun

 zeroth <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     delta = -1)$zeroth, theta)
 all.equal(zeroth, cumfun(theta))

 first <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     delta = -1, deriv = 1)$first, theta)
 all.equal(first, moofun(theta))

 second <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     delta = -1, deriv = 2)$second, theta)
 all.equal(second, voofun(theta))

 third <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     delta = -1, deriv = 3)$third, theta)
 all.equal(third, thoofun(theta))

 foo <- link(1, fam.zero.truncated.poisson(), delta = -1, deriv = 1)
 # does not matter what zeroth is, theta is not identifiable for
 # degenerate model
 is.finite(foo$zeroth)
 all.equal(0, foo$first)

 ##### multinomial #####

 set.seed(42)

 d <- 4

 theta <- matrix(rnorm(d * 25), ncol = d)

 cumfun <- function(theta, delta) {
     stopifnot(is.numeric(theta))
     stopifnot(is.finite(theta))
     stopifnot(length(theta) == d)
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == d)
     stopifnot(delta <= 0)
     stopifnot(any(delta == 0))
     inies <- delta == 0
     d.too <- sum(inies)
     if (d.too == 1) return(theta[inies])
     theta.too <- theta[inies]
     return(cumulant(theta.too, fam.multinomial(d.too))$zeroth)
 }

 moofun <- function(theta, delta) {
     stopifnot(is.numeric(theta))
     stopifnot(is.finite(theta))
     stopifnot(length(theta) == d)
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == d)
     stopifnot(delta <= 0)
     stopifnot(any(delta == 0))
     inies <- delta == 0
     d.too <- sum(inies)
     if (d.too == 1) as.numeric(inies)
     theta.too <- theta[inies]
     foo <- cumulant(theta.too, fam.multinomial(d.too), deriv = 1)$first
     bar <- rep(0, d)
     bar[inies] <- foo
     return(bar)
 }

 voofun <- function(theta, delta) {
     stopifnot(is.numeric(theta))
     stopifnot(is.finite(theta))
     stopifnot(length(theta) == d)
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == d)
     stopifnot(delta <= 0)
     stopifnot(any(delta == 0))
     inies <- delta == 0
     d.too <- sum(inies)
     if (d.too == 1) return(matrix(0, d, d))
     theta.too <- theta[inies]
     foo <- cumulant(theta.too, fam.multinomial(d.too), deriv = 2)$second
     bar <- matrix(0, d, d)
     baz <- matrix(0, d, d.too)
     baz[inies, ] <- foo
     bar[ , inies] <- baz
     return(bar)
 }

 thoofun <- function(theta, delta) {
     stopifnot(is.numeric(theta))
     stopifnot(is.finite(theta))
     stopifnot(length(theta) == d)
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == d)
     stopifnot(delta <= 0)
     stopifnot(any(delta == 0))
     inies <- delta == 0
     d.too <- sum(inies)
     if (d.too == 1) return(array(0, rep(d, 3)))
     theta.too <- theta[inies]
     foo <- cumulant(theta.too, fam.multinomial(d.too), deriv = 3)$third
     bar <- array(0, rep(d, 3))
     baz <- array(0, c(d, d, d.too))
     qux <- array(0, c(d, d.too, d.too))
     qux[inies, , ] <- foo
     baz[ , inies, ] <- qux
     bar[ , , inies] <- baz
     return(bar)
 }

 ##### one extra degeneracy #####

 delta <- c(0, 0, 0, -1)

 zeroth <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta)$zeroth)
 my.zeroth <- apply(theta, 1, cumfun, delta = delta)
 all.equal(zeroth, my.zeroth)

 first <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 my.first <- apply(theta, 1, moofun, delta = delta)
 all.equal(first, my.first)

 second <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 2)$second)
 my.second <- apply(theta, 1, voofun, delta = delta)
 all.equal(second, my.second)

 third <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 3)$third)
 my.third <- apply(theta, 1, thoofun, delta = delta)
 all.equal(third, my.third)

 link(first[ , 1], fam.multinomial(d), delta = delta, deriv = 1)

 ##### two extra degeneracies #####

 delta <- c(0, 0, -1, -2)

 zeroth <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta)$zeroth)
 my.zeroth <- apply(theta, 1, cumfun, delta = delta)
 all.equal(zeroth, my.zeroth)

 first <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 my.first <- apply(theta, 1, moofun, delta = delta)
 all.equal(first, my.first)

 second <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 2)$second)
 my.second <- apply(theta, 1, voofun, delta = delta)
 all.equal(second, my.second)

 third <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 3)$third)
 my.third <- apply(theta, 1, thoofun, delta = delta)
 all.equal(third, my.third)

 link(first[ , 1], fam.multinomial(d), delta = delta, deriv = 1)

 ##### fully degenerate #####

 delta <- c(0, -0.5, -1, -2)

 zeroth <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta)$zeroth)
 my.zeroth <- apply(theta, 1, cumfun, delta = delta)
 all.equal(zeroth, my.zeroth)

 first <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 my.first <- apply(theta, 1, moofun, delta = delta)
 all.equal(first, my.first)

 second <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 2)$second)
 my.second <- apply(theta, 1, voofun, delta = delta)
 all.equal(second, my.second)

 third <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 3)$third)
 my.third <- apply(theta, 1, thoofun, delta = delta)
 all.equal(third, my.third)

 ##### now link #####

 linkfun <- function(xi, delta) {
     stopifnot(is.numeric(xi))
     stopifnot(is.finite(xi))
     stopifnot(length(xi) == d)
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == d)
     stopifnot(delta <= 0)
     stopifnot(any(delta == 0))
     inies <- delta == 0
     d.too <- sum(inies)
     if (d.too == 1) return(rep(0, d))
     xi.too <- xi[inies]
     foo <- link(xi.too, fam.multinomial(d.too))$zeroth
     bar <- rep(0, d)
     bar[inies] <- foo
     return(bar)
 }

 dlinkfun <- function(xi, delta) {
     stopifnot(is.numeric(xi))
     stopifnot(is.finite(xi))
     stopifnot(length(xi) == d)
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == d)
     stopifnot(delta <= 0)
     stopifnot(any(delta == 0))
     inies <- delta == 0
     d.too <- sum(inies)
     if (d.too == 1) return(matrix(0, d, d))
     xi.too <- xi[inies]
     foo <- link(xi.too, fam.multinomial(d.too), deriv = 1)$first
     bar <- matrix(0, d, d)
     baz <- matrix(0, d, d.too)
     baz[inies, ] <- foo
     bar[ , inies] <- baz
     return(bar)
     bar[inies] <- foo
     return(bar)
 }

 ##### one extra degeneracy #####

 delta <- c(0, 0, 0, -1)

 xi <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 xi <- t(xi)

 zeroth <- apply(xi, 1, function(xi) link(xi, fam.multinomial(d),
     delta = delta)$zeroth)
 my.zeroth <- apply(xi, 1, linkfun, delta = delta)
 all.equal(zeroth, my.zeroth)

 first <- apply(xi, 1, function(xi) link(xi, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 my.first <- apply(xi, 1, dlinkfun, delta = delta)
 all.equal(first, my.first)

 ##### two extra degeneracies #####

 delta <- c(0, 0, -1, -2)

 xi <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 xi <- t(xi)

 zeroth <- apply(xi, 1, function(xi) link(xi, fam.multinomial(d),
     delta = delta)$zeroth)
 my.zeroth <- apply(xi, 1, linkfun, delta = delta)
 all.equal(zeroth, my.zeroth)

 first <- apply(xi, 1, function(xi) link(xi, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 my.first <- apply(xi, 1, dlinkfun, delta = delta)
 all.equal(first, my.first)

 ##### fully degenerate #####

 delta <- c(0, -0.5, -1, -2)

 xi <- apply(theta, 1, function(theta) cumulant(theta, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 xi <- t(xi)

 zeroth <- apply(xi, 1, function(xi) link(xi, fam.multinomial(d),
     delta = delta)$zeroth)
 my.zeroth <- apply(xi, 1, linkfun, delta = delta)
 all.equal(zeroth, my.zeroth)

 first <- apply(xi, 1, function(xi) link(xi, fam.multinomial(d),
     delta = delta, deriv = 1)$first)
 my.first <- apply(xi, 1, dlinkfun, delta = delta)
 all.equal(first, my.first)

