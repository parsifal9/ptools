# interpolate survival function output by cross validation
# assumes no ties in survival times so can use a common set of times
# i.e for each fold ties are not broken in a different way

survinterp <- function(stime, xsurv)
{
	stime <- sort(stime)
	n <- ncol(xsurv)

	for (i in 1:n) {
		tt <- is.na(xsurv[, i])
		res <- approx(stime[!tt], xsurv[!tt, i], xout = stime[tt], rule = 2)
		xsurv[tt, i] <- res$y
	}

	list(stime = stime, xsurv = xsurv)
}
