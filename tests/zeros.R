
 library(aster2)

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

 czero <- function(data) {
     validasterdata(data)
     nnode <- length(data)
     .C(aster2:::C_aster_predecessor_zero_cond,
         nnode = as.integer(nnode), pred = as.integer(data$repred),
         resp = as.double(data$redata[[data$response.name]]),
         result = logical(nnode))$result
 }

 uzero <- function(data) {
     validasterdata(data)
     nnode <- length(data)
     aster2:::fam.clear()
     for (i in seq(along = data$families))
         aster2:::fam.set(data$families[[i]])
     out <- .C(aster2:::C_aster_predecessor_zero_unco,
         nnode = as.integer(nnode), pred = as.integer(data$repred),
         group = as.integer(data$regroup), code = as.integer(data$recode),
         delta = as.double(data$redelta), result = logical(nnode))
     aster2:::fam.clear()
     return(out$result)
 }

 sally <- czero(fred)

 resp <- fred$redata[[fred$response.name]]
 pred <- fred$repred
 ypred <- c(NA, resp)[pred + 1]
 ypred[is.na(ypred)] <- fred$initial[is.na(ypred)]
 
 identical(sally, ypred == 0)

 sally <- uzero(fred)
 all(sally == FALSE)

 varb <- as.character(fred$redata[[fred$varb.name]])
 unique(varb)

 set.seed(43)

 delta <- fred$redelta
 all(delta == 0)
 delta[grepl("[mbp]", varb) & resp == 0 & runif(length(fred)) < 1 / 5] <- -1
 fred$redelta <- delta
 is.validasterdata(fred)
 sum(delta != 0)
 
 sally <- czero(fred)
 identical(sally, ypred == 0)

 sally <- uzero(fred)
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
 sum(my.sally == 2)
 sum(my.sally == 1)
 my.sally <- my.sally == 1
 identical(sally, my.sally)

