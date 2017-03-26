
 library(aster2)

 set.seed(42)

 data(test1)

 fred <- asterdata(test1,
     vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
     pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
     code = c(1, 1, 1, 2, 2, 3, 4, 5),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))
 is.validasterdata(fred)
 names(fred)
 names(fred$redata)

 theta <- rnorm(length(fred))
 try(validtheta(fred, theta))
 try(validtheta(fred, theta, mod = "cond"))

 varb <- as.character(fred$redata[[fred$varb.name]]) 
 pred <- fred$repred
 resp <- fred$redata[[fred$response.name]]
 init <- fred$initial
 ypred <- ifelse(pred == 0, init, c(NA, resp)[pred + 1])

 theta.cond <- theta
 theta.unco <- theta
 theta.unco[varb == "n2"] <- (- abs(theta.unco[varb == "n2"]))
 theta.cond[varb == "n2" & ypred > 0] <-
     (- abs(theta.unco[varb == "n2" & ypred > 0]))

 is.validtheta(fred, theta.unco)
 is.validtheta(fred, theta.cond, mod = "cond")
 
 xi <- rnorm(length(fred))
 try(validxi(fred, xi))
 try(validxi(fred, xi, mod = "cond"))

 fixup <- function(xi, varb, ypred) {
     stopifnot(length(xi) == length(varb))
     stopifnot(is.numeric(xi))
     stopifnot(is.finite(xi))
     stopifnot(is.character(varb))
     if (missing(ypred)) {
         ypred <- rep(1, length(xi))
     } else {
         stopifnot(is.numeric(ypred))
         stopifnot(length(ypred) == length(xi))
     }
     todo <- ypred != 0
     xi[todo & varb == "z1"] <- 1 + abs(xi[todo & varb == "z1"])
     xi[todo & varb == "p1"] <- abs(xi[todo & varb == "p1"])
     xi[todo & varb == "b1"] <- 1 / (1 + abs(xi[todo & varb == "b1"]))
     xi[todo & varb == "n2"] <- 1 + xi[todo & varb == "n1"]^2
     my.xi.m1 <- exp(xi[todo & varb == "m1"])
     my.xi.m2 <- exp(xi[todo & varb == "m2"])
     my.xi.m3 <- exp(xi[todo & varb == "m3"])
     xi[todo & varb == "m1"] <- my.xi.m1 / (my.xi.m1 + my.xi.m2 + my.xi.m3)
     xi[todo & varb == "m2"] <- my.xi.m2 / (my.xi.m1 + my.xi.m2 + my.xi.m3)
     xi[todo & varb == "m3"] <- my.xi.m3 / (my.xi.m1 + my.xi.m2 + my.xi.m3)
     return(xi)
 }

 xi.unco <- fixup(xi, varb)
 xi.cond <- fixup(xi, varb, ypred)

 is.validxi(fred, xi.unco)
 is.validxi(fred, xi.cond, mod = "cond")

 ########## now with direction of recession ##########

 delta <- fred$redelta
 all(delta == 0)
 delta[grepl("[mbp]", varb) & resp == 0 & runif(length(fred)) < 1 / 5] <- -1
 fred$redelta <- delta
 is.validasterdata(fred)
 sum(delta != 0)
 
 try(validtheta(fred, theta))
 try(validtheta(fred, theta, mod = "cond"))

 is.validtheta(fred, theta.unco)
 is.validtheta(fred, theta.cond, mod = "cond")
 
 my.sally <- 2 * as.numeric(delta < 0)
 # when loop done
 # my.sally == 2 means delta < 0 and predecessor NOT a. s. zero
 # my.sally == 1 means predecessor a. s. zero (regardless of delta)
 repeat {
     save.my.sally <- my.sally
     foom <- c(0, my.sally)[pred + 1]
     my.sally[foom != 0] <- 1
     if (identical(my.sally, save.my.sally)) break
 }
 my.sally <- my.sally == 1

 try(validxi(fred, xi))
 try(validxi(fred, xi, mod = "cond"))

 fixup <- function(xi, varb, ypred, delta) {
     stopifnot(is.character(varb))
     stopifnot(length(xi) == length(varb))
     stopifnot(is.numeric(xi))
     stopifnot(is.finite(xi))
     stopifnot(is.numeric(ypred))
     stopifnot(is.finite(ypred))
     stopifnot(length(ypred) == length(varb))
     stopifnot(is.numeric(delta))
     stopifnot(is.finite(delta))
     stopifnot(length(delta) == length(varb))
     todo <- ypred != 0
     lower <- delta < 0
     xi[todo & varb == "z1"] <- 1 + abs(xi[todo & varb == "z1"])
     xi[todo & varb == "z1" & lower] <- 1
     xi[todo & varb == "p1"] <- abs(xi[todo & varb == "p1"])
     xi[todo & varb == "p1" & lower] <- 0
     xi[todo & varb == "b1"] <- 1 / (1 + abs(xi[todo & varb == "b1"]))
     xi[todo & varb == "b1" & lower] <- 0
     xi[todo & varb == "n2"] <- 1 + xi[todo & varb == "n1"]^2
     my.xi <- ifelse(lower, 0, exp(xi))
     my.xi.m1 <- my.xi[todo & varb == "m1"]
     my.xi.m2 <- my.xi[todo & varb == "m2"]
     my.xi.m3 <- my.xi[todo & varb == "m3"]
     xi[todo & varb == "m1"] <- my.xi.m1 / (my.xi.m1 + my.xi.m2 + my.xi.m3)
     xi[todo & varb == "m2"] <- my.xi.m2 / (my.xi.m1 + my.xi.m2 + my.xi.m3)
     xi[todo & varb == "m3"] <- my.xi.m3 / (my.xi.m1 + my.xi.m2 + my.xi.m3)
     return(xi)
 }

 xi.unco <- fixup(xi, varb, as.numeric(! my.sally), delta)
 xi.cond <- fixup(xi, varb, ypred, delta)

 is.validxi(fred, xi.unco)
 is.validxi(fred, xi.cond, mod = "cond")

