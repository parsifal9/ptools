\name{xvalidate}

\alias{xvalidate}

\title{Produces cross validated fitted values}

\description{
Produces cross validated fitted values for methods \code{\link{sdda}},
\code{\link{HGmultc}}, \code{\link{HGordcat}}, \code{\link{HGgaussian}},
\code{\link{HGglm}}, \code{\link{HGcoxreg}}, \code{\link{HGsurv}}.
}

\usage{
xvalidate(X, y, event = NULL, method, fold = NULL, trace = FALSE,
    weights = rep(1, nrow(X)), \dots)
}

\arguments{

\item{X}{real matrix with n rows (samples) and p columns (variables or genes)}

\item{y}{
	real n by 1 vector of response values appropriate to the method. For
	\code{\link{HGcoxreg}} and \code{\link{HGsurv}}, \code{y} is the vector of
	censored survival times
}

\item{event}{
	Only required for methods \code{\link{HGcoxreg}} and \code{\link{HGsurv}}.
	Real n by 1 censoring indicator - 1 denotes survival time uncensored,
	0 - censored
}

\item{method}{
	one of \code{\link{sdda}}, \code{\link{HGmultc}}, \code{\link{HGordcat}},
	\code{\link{HGgaussian}}, \code{\link{HGglm}}, \code{\link{HGcoxreg}},
	\code{\link{HGsurv}}
}

\item{fold}{
	an integer specifying the number of folds to use in the cross validation.
	Default value NULL gives leave one out cross validation
}

\item{trace}{
	a logical variable with value TRUE if information about cross validation
	progress to be output and value FALSE if no output required.
}

\item{weights}{real n by 1 vector of observation weights, default to 1}

\item{\dots}{
	other parameters appropriate to the selected method. Not specifying these
	results in the default values for these parameters
}

}

\details{
The data is divided into fold groups. The fitted value for each group is
obtained by using a model built from the data for the remaining groups.

The parameters \code{bbess}, \code{kbess} and the \code{weights} to be used
during the fitting can be specified. They default to give the Jeffreys
hyperprior and unit weights respectively.

For all methods except \code{\link{sdda}}, \code{\link{HGmultc}} and
\code{\link{HGordcat}}, \code{b0sc} can be specified.

In adddition, for the methods \code{\link{HGgaussian}} and
\code{\link{HGglm}}, the \code{scale} parameter can be fixed or estimated.

Except for \code{\link{sdda}}, \code{\link{HGmultc}} and
\code{\link{HGordcat}}, the cross validated fitted value is defined as the
models fitted value. This may not always be what you want, e.g. for logistic
regression using a glm.

For \code{\link{HGmultc}}, \code{\link{HGordcat}} and \code{\link{sdda}} the
usual cross validated class labels are produced. For \code{\link{sdda}},
\code{weights} are ignored.

For \code{\link{HGsurv}} and \code{\link{HGcoxreg}} the fitted
values are defined to be the linear predictor values.
}

\value{ An n by 1 vector of cross validated fitted values, except for methods HGcoxreg and HGsurv
in which case it is a list with components
  \item{xval.lp }{cross validated linear predictor}
  \item{xval.surv}{an n by n+2 matrix of cross validated survival curves, with rows indexed by time. The first two
  columns correspond to the sorted survival time vector and the correspondingly sorted event vector.}
}

\author{ Glenn Stone and Harri Kiiveri}

\note{ Missing values in \code{X} and \code{y} are not allowed }

\section{Warning }{If a method is not specified an error occurs. There is no type checking of arguments. }

\seealso{
\code{\link{xval.survival}},
\code{\link{HGglm}},
\code{\link{HGordcat}},
\code{\link{HGmultc}},
\code{\link{HGgaussian}},
\code{\link{HGcoxreg}},
\code{\link{HGsurv}},
\code{\link{DeleteRepeat}}
}

\examples{
# EXAMPLE - generate classification data and run xvalidate

x <- matrix(rnorm(200 * 200), ncol = 200)
y <- (x[, 10] > 0) + 1

# 10 fold cross validation
res <- xvalidate(x, y, method = HGmultc, fold = 10)

# plot y against cross validated fitted value
plot(y, res)
}

\keyword{models}
\keyword{multivariate}

