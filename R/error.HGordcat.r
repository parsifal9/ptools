error.HGordcat<- function (object, X, y, weights, fv, ...)
{
	if (missing(fv)) {
		#plugin case
		b <- tapply(weights, list(y, predict.HGordcat(object, X, type = "class")),
		simplify = TRUE, FUN = "sum")
		e <- 1 - sum(diag(b)) / sum(weights)
	}
	else {
		# cross validation case
		e <- sum(weights * (y != fv)) / sum(weights)
	}
	e
}
