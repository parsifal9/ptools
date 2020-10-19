HGcontrol <- function(tolerance = 1e-4, epsilon = 1e-4, maxit = 200,
	tolc = 1e-2, eta.lim = 100)
{
# tolerance - convergence parameter
# epsilon - parameter for thresholding estimates
# maxit - maximum number of E-steps
# tolc - line search parameter
# eta.lim - maximum absolute value of linear predictor before truncation

	if (!is.numeric(tolerance) || tolerance <= 0)
		stop("value of 'tolerance' must be > 0")
	if (!is.numeric(epsilon) || epsilon <= 0)
		stop("value of 'epsilon' must be > 0")
	if (!is.numeric(maxit) || maxit <= 0)
		stop("maximum number of iterations must be > 0")
	if (!is.numeric(tolc) || tolc <= 0)
		stop("value of 'tolc' must be > 0")
  if (!is.numeric(eta.lim) || eta.lim <= 0)
		stop("value of 'eta.lim' must be > 0")

	list(tolerance = tolerance, epsilon = epsilon, maxit = maxit, tolc = tolc,
       eta.lim=eta.lim)
}
