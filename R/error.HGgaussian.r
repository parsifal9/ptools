error.HGgaussian <- function (object, X, y, weights, fv, ...)
{
	if (missing(fv))
		r <- (y - object$fv)
	else
		r <- (y - fv)

	mean(r * r * weights)
}
