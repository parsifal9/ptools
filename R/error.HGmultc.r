error.HGmultc <- function (object, X, y, weights, fv, ...)
{
    if (missing(weights)) {
        weights <- rep(1, dim(X)[[2]])
    }
    if (missing(fv)) {
       fv <- predict.HGmultc(object,X, type = "class")  # define fv in a separate line
       e <- sum(weights * (y != fv))/sum(weights)
    }
    else {
        e <- sum(weights * (y != fv))/sum(weights)
    }
   e
}
