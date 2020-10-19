error.HGsvm <- function (object, X, y, weights, fv, ...)
{
	if (missing(weights)) {
		weights <- rep(1, dim(X)[[2]])
	}

	if (missing(fv)) {
		b <- tapply(weights, list(y, predict.HGsvm(object, X)),
			simplify = TRUE, FUN = "sum")
		e <- 1 - sum(diag(b)) / sum(weights)
	}
	else {
		e <- sum(weights * (y != fv)) / sum(weights)
	}
	e
}
