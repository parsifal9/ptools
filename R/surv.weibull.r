#
# functions to implement weibull in HGsurv
#

loginthaz <- function (logt, alpha)
{
	alpha * logt
}

dloginthaz <- function (logt, alpha)
{
	alpha / exp(logt)
}

aupdate <- function (f, logt, alpha0, p0, weights)
{
	alpha1 <- sum(logt * (p0 - f) * weights)
	alpha1 <- sum(weights * f) / alpha1
	alpha0 + 0.1 * (alpha1 - alpha0)
}

ainit <- function ()
{
	aa <- 1
	aa
}
