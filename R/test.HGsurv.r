test.HGsurv <- function(testing = TRUE, test.name = "HGsurv")
{
	set.seed(1)

	m   <- 20
	n   <- 35
	x   <- matrix(rnorm(m * n), nrow = m)  # Generate data
	lp  <- 3 + x[, 1]
	a   <- 0.5                         # Weibull shape parameter 0.5
	t   <- rweibull(m, a, exp(-lp/a))  # Simulate survival times
	c   <- rep(0:1, length.out = m)    # No censoring
	try <- HGsurv(x, t, c, b0sc = 0.5) # Fit model

	test.compare(try, truth, test.name, testing)
}
