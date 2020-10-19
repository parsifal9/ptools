\name{gethyperpars}

\alias{gethyperpars}

\title{
Determines a table of crossvalidated error rates indexed by values of the
hyperparameters bbess, kbess and b0sc
}

\description{
Determines a table of crossvalidated error rates indexed by values of the
hyperparameters \code{bbess}, \code{kbess} and \code{b0sc}.
}

\usage{
gethyperpars(X, y, event = NULL, method, xvalid = TRUE,
    fold = NULL, trace = FALSE, weights = rep(1, nrow(X)),
    bbess = 1e7, kbess = 0, b0sc = c(1, 5, 10, 15, 20, 25), \dots)
}

\arguments{
\item{X}{
	real matrix with n rows (samples) and p columns (variables or genes)
}

\item{y}{
	real n by 1 vector of response values appropriate to the method. For
	\code{\link{HGcoxreg}} and \code{\link{HGsurv}}, \code{y} is the vector
	of censored survival times
}

\item{event}{
	Only required for \code{\link{HGcoxreg}} and \code{\link{HGsurv}}.
	Real n by 1 censoring indicator - 1 denotes survival time uncensored,
	0 - censored
}

\item{method}{
	one of \code{\link{HGmultc}}, \code{\link{HGordcat}},
	\code{\link{HGgaussian}}, \code{\link{HGglm}}, \code{\link{HGcoxreg}},
	\code{\link{HGsurv}}
}

\item{xvalid}{
	a logical variable with value TRUE if cross validation is required,
	FALSE	if plugin error estimates are required. Default is TRUE.
}

\item{fold}{
	an integer specifying the number of folds to use in the cross validation.
	Default value NULL gives leave one out cross validation.
}

\item{trace}{
	a logical variable with value TRUE if information about cross validation
	progress to be output and value FALSE if no output required.
}

\item{weights}{real n by 1 vector of observation weights, default to 1.}

\item{bbess}{
	vector of test values of scale parameter in prior for HG method,
	default 1e7.
}

\item{kbess}{
	vector of test values of shape parameter in prior for HG method - 0 to 1,
	default 0
}

\item{b0sc}{
	vector of test values of length of initial beta vector - likely useful
	values in the range 1 - 25. Ignored in the methods \code{\link{HGmultc}},
	\code{\link{HGordcat}}.
}

\item{\dots}{ any other options appropriate to the chosen method. }
}

\details{
Generates (crossvalidated if specified) error rates for all combinations of
\code{bbess}, \code{kbess} and \code{b0sc} specified in the input arguments.
Results are returned in a (two) three dimensional array. The default values
of \code{bbess} and \code{kbess} give the Jeffreys hyperprior.

For all methods except \code{\link{HGmultc}} and \code{\link{HGordcat}},
a vector of values of \code{b0sc} can be specified. They are ignored for
\code{\link{HGmultc}} and \code{\link{HGordcat}} and a two dimensional array
is produced.

In adddition, for the methods \code{\link{HGgaussian}} and \code{HGglm},
the \code{scale} parameter can be fixed or estimated.

Except for \code{\link{HGmultc}} and \code{\link{HGordcat}}, the cross
validated error is defined as the weighted mean square error calculated
from \code{y} and the models fitted value. This may not always be what you
want, e.g. for logistic regression using a glm.

For \code{\link{HGmultc}} and \code{\link{HGordcat}} the usual weighted
misclassification error rate is produced. NOTE: For \code{\link{HGsurv}}
and \code{\link{HGcoxreg}} error rates are currently returned as NA.
}

\value{ A list with components

\item{xvtable}{
	A (two) three dimensional array of error rates  of dimension (p by) q by r,
	where p is the length of the \code{b0sc} vector, q is the length of the
	\code{bbess} vector and r is the length of the \code{kbess} vector.
}
}

\author{Harri Kiiveri}

\note{ Missing values in \code{X} and \code{y} are not allowed }

\section{Warning}{
If a method is not specified an error occurs. There is no type checking of
arguments.
}

\seealso{
\code{\link{HGglm}},
\code{\link{HGordcat}},
\code{\link{HGmultc}},
\code{\link{HGgaussian}},
\code{\link{HGcoxreg}},
\code{\link{HGsurv}},
\code{\link{DeleteRepeat}}
}

\examples{
# EXAMPLE - generate linear regression data and run gethyperpars

x <- matrix(rnorm(200 * 200), ncol = 200)
y <- (x[, 10] > 0) + 1

res <- gethyperpars(x, y, method = HGmultc, xvalid = TRUE, fold = 10,
    kbess = c(0, 0.5, 1), b0sc = c(1, 5, 10))

# Look at table
res$xvtable
}

\keyword{ models }
\keyword{multivariate}
