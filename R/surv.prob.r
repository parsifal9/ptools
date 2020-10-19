#  get survival probabilites for a vector t0 of specified times
#  one for each individual
# t0 can be a scalar
surv.prob <- function(stime, xsurv, t0)
{
# computes probability that survival times exceed t0 using
# the individual survival curves in the cols of xsurv
	n <- ncol(xsurv)

	if (length(t0) == 1)
		t0 <- rep(t0, n)

	stime <- sort(stime)
	phat <- rep(0, n)

	for (i in 1:n) {
		res <- approx(stime, xsurv[, i], xout = t0[i], rule = 2)
		phat[i] <- res$y
	}

	phat
}
