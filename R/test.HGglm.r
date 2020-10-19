test.HGglm <- function(testing = TRUE, test.name = "HGglm")
{
	set.seed(1)

##
## EXAMPLE 1 - Logistic regression using HGglm
##

	n <- 50
	w <- abs(sin(1:n))
	x <- matrix(rnorm(n * n), nrow = n)
	y <- as.numeric(x[, 1] > 0)

	res1a <- HGglm(x, y)              # Fit GLM with default weights and default Jeffreys hyperprior
	res1b <- HGglm(x, y, kbess = 1)   # Use Lasso prior with default scale parameters
	res1c <- HGglm(x, y, weights = w) # Use weights in w

	## Check results
	try <- list(res1a = res1a, res1b = res1b, res1c = res1c)

	test.compare(try, truth, test.name, testing)
}
