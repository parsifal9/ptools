error.HGL1reg <- function (object, X, y, weights, fv, ...)
{
	if (missing(fv))
		r <- (y - object$fv)
	else
		r <- (y - fv)

	mean(abs( r) * weights)
}
