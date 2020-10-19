error.HGglm <- function (object, X, y, weights, fv, ...)
{
	if (object$model == "B")
		y <- y / object$bd

	if (missing(fv))
		r <- (y - object$fv)
	else
		r <- (y - fv)

	mean(r * r * weights)
}
