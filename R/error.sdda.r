error.sdda <- function (object, X, y, weights, fv, ...)
{
# weights ignored here
	if (missing(fv)) {
		#plugin case
		tt <- sum(diag(table(y, predict(object, X, type = "class"))))
		e <- 1 - as.double(tt) / length(y)
	}
	else {
		# cross validation case
		e <- sum((y != fv)) / length(y)
	}

	e
}
