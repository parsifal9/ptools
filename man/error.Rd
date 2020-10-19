\name{error}
\alias{error}
\alias{error.fastRP}
\alias{error.HGgaussian}
\alias{error.HGglm}
\alias{error.HGmultc}
\alias{error.HGsurv}
\alias{error.sdda}

\title{Produces an error estimate}

\description{
Produces an error estimate. Methods exist for the class \code{\link{fastRP}},
\code{\link{HGgaussian}}, \code{\link{HGglm}}, \code{\link{HGmultc}},
\code{\link{HGsurv}} and \code{\link{sdda}}.
}

\usage{
error(object, X, y, weights, fv, \dots)
\method{error}{fastRP}(object, X, y, weights, fv, \dots)
\method{error}{HGgaussian}(object, X, y, weights, fv, \dots)
\method{error}{HGglm}(object, X, y, weights, fv, \dots)
\method{error}{HGmultc}(object, X, y, weights, fv, \dots)
\method{error}{HGsurv}(object, X, y, weights, fv, \dots)
\method{error}{sdda}(object, X, y, weights, fv, \dots)
}

\arguments{
\item{object}{fitted model}
\item{X}{
	real matrix with n rows (samples) and p columns (variables or genes).
	The matrix should have the correct number of columns
}
\item{y}{
	real n by 1 vector of response values appropriate to the method. For
	\code{\link{HGcoxreg}} and \code{\link{HGsurv}}, \code{y} is the vector
	of censored survival times
}
\item{weights}{real n by 1 vector of observation weights, default to 1}
\item{fv}{
	The fitted values from the model (object) applied to the data (\code{X}
	and \code{y})
}
\item{\dots}{Any other arguments, not currently used}
}

\value{
}

\details{
The function calcuates the error rate. \code{X} and \code{y} can be the data
the model was fitted on or they can be new data (of the correct size)
}

\examples{
set.seed(123)
w <- abs(sin(1:200))
x <- matrix(rnorm(200 * 200), nrow = 200)
y <- (x[, 1] > 0) + 1

# use default weights and default Jeffreys hyperprior
res <- HGmultc(x, y)

# display misclassification table
table(y, res$class)

error(res, x, y)

# make the problem a bit harder

set.seed(123)
y[31:50] <- sample(c(1, 2), 20, replace = TRUE)

res <- HGmultc(x, y)
table(y, res$class)
error(res, x, y)

# get the error rate on new data
set.seed(123)
newdata <- matrix(rnorm(200 * 200), nrow = 200)
newy <- (newdata[, 1] > 0) + 1

error(res, newdata, newy)
}
\author{Rob Dunne}
\keyword{}
