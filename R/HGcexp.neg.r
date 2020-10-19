HGcexp.neg <- function(bbess, kbess, beta)
{
	if (length(beta) == 0)
		return()

	beta <- matrix(beta, ncol = 1)
	# computes expected value of 1/sigma2 for the NEG prior
	lambda <- kbess

	if (lambda > 1000) {
		warning("lambda exceeds 1000")
	}

	g <- 1 / bbess
	beta <- abs(beta) / g
	n <- length(beta)
	val <- rep(0, n)

	vn <- -2 * lambda - 2
	vd <- vn + 1
	const <- -vd / g^2

	val <- apply(beta, 1, FUN = expval, vn, vd, const)
	val
}

expval <- function(beta, vn, vd, const)
{
	#dn <- paracylfunc(vn,beta)
	#dd <- paracylfunc(vd,beta)

	#use asymptotic approx for ratio if numbers underflow
	if (abs(beta) > 40) {
		r <- 1 / beta
	}
	else {
		dn <- paracylfunc(vn, beta)
		dd <- paracylfunc(vd, beta)

		if (is.nan(dd$pdf) | (dd$pdf == 0)) {
			r <- 1 / beta
		}
		else{
			r <- dn$pdf / dd$pdf
		}
	}

	val <- const * r / beta
	val
}
