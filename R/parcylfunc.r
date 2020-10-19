paracylfunc <- function(v = 5.5, x = 10, pdf = 0, pdd = 0)
{
	# D(v,x) is in pdf, its derivative is in pdd
	nv <- abs(trunc(v))
	nv <- max(nv + 1, 1)
	dp <- rep(0, nv)
	dv <- dp

	res <- .Fortran("pbdv", v = v, x = x, dv = dv, dp = dp, pdf = pdf, pdd = pdd,
		PACKAGE = "RChip")

	list(v = res$v, x = res$x, dv = res$dv, dp = res$dp, pdf = res$pdf, pdd = res$pdd)
}
