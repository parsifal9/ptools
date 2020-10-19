\name{HGcontrol}
\alias{HGcontrol}

\title{ Auxiliary for Controlling HG Fitting }
\description{
Auxiliary function as user interface for HG fitting. Typically only
used when calling the HG* functions.
}

\usage{
HGcontrol(tolerance = 1e-04, epsilon = 1e-04, maxit = 200, tolc = 1e-2, eta.lim=100)
}

\arguments{
  \item{tolerance}{ convergence parameter. }
  \item{epsilon}{ parameter for thresholding estimates. }
  \item{maxit}{ integer giving the maximum number of E-steps. }
  \item{tolc}{ line search parameter}
  \item{eta.lim}{ maximum absolute value of linear predictor before truncation}
}

\details{
This function is based on \code{\link{glm.control}}.
}

\value{
A list with the arguments as components.
}

\author{ Tim Keighley }

\seealso{ \code{\link{HGengine}} }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(tolerance = 1e-4, epsilon = 1e-4, maxit = 200, tolc = 1e-2, eta.lim=100) {
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
if (!is.numeric(eta.lim) || eta.lim <= 0)
stop("eta.lim must be > 0")

list(tolerance = tolerance, epsilon = epsilon, maxit = maxit, eta.lim=eta.lim)
  }
}

\keyword{ misc }
