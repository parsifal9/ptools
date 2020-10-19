test.HGmultc <- function(testing = TRUE, test.name = "HGmultc")
{
	set.seed(1)

	n <- 20
	w <- abs(sin(1:n))
	x <- matrix(rnorm(n * n), nrow = n)
	y <- (x[, 1] > 0) + 1

	res1a <- HGmultc(x, y)              # Use default weights and default Jeffreys hyperprior
	res1b <- HGmultc(x, y, kbess = 1)   # Use Lasso prior with default scale parameter
	res1c <- HGmultc(x, y, weights = w) # Use weights in w

	## Check results
	try <- list(res1a = res1a, res1b = res1b, res1c = res1c)

	test.compare(try, truth, test.name, testing)
}
