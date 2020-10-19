test.HGgaussian <- function(testing = TRUE, test.name = "HGgaussian")
{
	set.seed(1)

	## EXAMPLE 1 - Generate regression data and test
	x <- matrix(rnorm(200 * 200), ncol = 200)
	y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)

# Scale fixed to 1
	res1a <- HGgaussian(x, y)
# Also estimate scale parameter
	res1b <- HGgaussian(x, y, sc = 1, b0sc = 1)
# Restart iterations from previous values
	res1c <- HGgaussian(x, y, sc = 1, b0sc = 1, initb = res1b$beta)

	## EXAMPLE 2 - Regression spline example
	x  <- seq(-2, 2, .02) + 1e-8
	y  <- sin(3 * x) / (3 * x)
	yy <- y + 0.125 * rnorm(201)

	tmp <- outer(x, x, FUN = "-")
	tmp <- tmp^3 * as.numeric(tmp > 0)
	xx  <- cbind(rep(1,201), x, x^2, x^3, tmp)  # Compute regression spline basis

	res2 <- HGgaussian(xx, yy)

	## Check results
	try <- list(res1a = res1a, res1b = res1b, res1c = res1c, res2 = res2)

	test.compare(try, truth, test.name, testing)
}
