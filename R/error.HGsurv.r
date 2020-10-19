error.HGsurv <- function (object, X, y, weights, fv, ...)
{
	if (missing(fv))
		r <- (y - object$lp)  #use linear predictor
	else
		r <- (y - fv)

	# mean(r * r * weights)
	return(NA)  #until get useful critera
}
