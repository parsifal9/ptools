error.HGcoxreg <- function (object, X, y, weights, fv, ...)
{
	if (missing(fv))
		r <- (y - object$lp)    #for now use linear predictor
	else
		r <- (y - fv)

	#mean(r * r * weights)
	return(NA) # until appropriate criterion found
}
