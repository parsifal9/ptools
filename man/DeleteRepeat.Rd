\name{DeleteRepeat}
\alias{DeleteRepeat}

\title{
Fits a sequence of HG (GeneRave) or sdda models by successively removing
selected sets of variables from consideration
}

\description{
Sequentially applies the selected method, with each application considering
only variables which have not been selected previously. Calculates plugin and
crossvalidated fitted values and error rates for each chosen model.
}

\usage{
DeleteRepeat(X, y, event = NULL, method = HGmultc, rounds = 1,
    xvalid = TRUE, fold = NULL, trace = FALSE,
    weights = rep(1, nrow(X)), \dots)
}

\arguments{
	\item{X}{
		real matrix with n rows (samples) and p columns (variables or genes)
	}
	\item{y}{
		real n by 1 vector of response values appropriate to the method.
		For \code{\link{HGcoxreg}} and \code{\link{HGsurv}}
		\code{y} is the vector of censored survival times
	}
	\item{event}{
		Only required for methods \code{\link{HGcoxreg}} and
		\code{\link{HGsurv}}. Real n by 1 censoring indicator - 1
		denotes survival time uncensored, 0 - censored
	}
	\item{method}{
		one of \code{\link{sdda}}, \code{\link{HGmultc}},
		\code{\link{HGordcat}}, \code{\link{HGgaussian}}, \code{\link{HGglm}},
		\code{\link{HGcoxreg}}, \code{\link{HGsurv}}
	}
	\item{rounds}{
		an integer specifying how many times to sequentially apply the HG
		method. Default is 2.
	}
	\item{xvalid}{
		a logical variable with value TRUE if cross validation is required,
		FALSE if plugin error estimates are required. Default is TRUE.
	}
	\item{fold}{
		an integer specifying the number of folds to use in the cross
		validation. Default value NULL gives leave one out cross validation.
	}
	\item{trace}{
		a logical variable with value TRUE if information about cross
		validation progress to be output and value FALSE if no output required.
	}
	\item{weights}{
		real n by 1 vector of observation weights, default to 1.
	}
	\item{\dots}{
		other parameters appropriate to the selected method. Not specifying
		these results in the default values for these parameters.
	}
}

\details{
The HG method is applied \code{rounds} times using only the variables which
have not been selected in previous rounds.

The parameters \code{bbess}, \code{kbess} and the \code{weights} to be used
during the fitting can be specified. They default to give the Jeffreys
hyperprior and unit weights respectively.

For all methods except \code{\link{sdda}}, \code{\link{HGmultc}} and
\code{\link{HGordcat}}, \code{b0sc} can be specified.

In adddition, for the methods \code{\link{HGgaussian}} and
\code{\link{HGglm}}, the \code{scale} parameter can be fixed or estimated.

Except for \code{\link{sdda}}, \code{\link{HGmultc}} and
\code{\link{HGordcat}}, the cross validated error is defined as the weighted
mean square error calculated from \code{y} and the models fitted value. This
may not always be what you want, e.g. for logistic regression using a glm.

For \code{\link{HGmultc}} and \code{\link{HGordcat}}, the usual weighted
misclassification error rate is produced. For \code{\link{sdda}}, weights are
ignored and the usual misclassification error rate is returned. For
\code{\link{HGsurv}} and \code{\link{HGcoxreg}} error rates are returned as NA.
See \code{\link{eval.surv.fit}} for a method for computing observed and expected
counts in specified time intervals using the output of this function.
}

\value{
A list with components
\item{genes.when}{
	a vector of length p. Components corresponding to variables selected at
	round k for k = 1, \ldots, \code{rounds} have value k. The remaining values
	are -1.
}
\item{error.rates}{
	a vector of length \code{rounds} containing the error rate for the model
	chosen at each round.
}
\item{xvfv}{
	a real matrix of dimension n by \code{rounds} with k-th column the
	(crossvalidated if \code{xvalid=TRUE}) fitted values for the model
	selected during the k-th round. Currently, for the methods
	\code{\link{HGsurv}} and \code{\link{HGcoxreg}} the (cross validated if
	\code{xvalid=TRUE}) linear predictor is returned. For
	\code{\link{HGmultc}} and \code{\link{HGordcat}},
	(cross validated if \code{xvalid=TRUE}) fitted classes are returned. }
\item{grp}{Only if xvalid=TRUE. An n by rounds matrix with ith column
  identifying the fold each observation belonged to for the ith round.}
\item{xval.surv}{A NULL list unless method is HGcoxreg or HGsurv, in which
  cases it is a list with ith component being an n by n+2 matrix of fitted
  (cross validated,if xvalid=T) survival curves ,with rows indexed by time,
  for the model chosen in the ith round. The first two columns of these matrices
  are the sorted survival times and correspondingly sorted event vector.}
}

%\references{ }

\author{ Glenn Stone and Harri Kiiveri}

\note{ Missing values in \code{X} and \code{y} are not allowed }

\section{Warning}{
	If a method is not specified an error occurs. There is no type checking of
	arguments.
}

\seealso{
\code{\link{HGcoxreg}},
\code{\link{HGgaussian}},
\code{\link{HGglm}},
\code{\link{HGmultc}},
\code{\link{HGordcat}},
\code{\link{HGsurv}}
}

\examples{
# EXAMPLE - generate linear regression data and run DeleteRepeat

x <- matrix(rnorm(200 * 200), ncol = 200)
y <- (x[, 10] > 0) + 1

# 2 rounds (default) no cross validation
res <- DeleteRepeat(x, y, method = HGmultc, xvalid = FALSE, rounds = 2)

# variables chosen on first round
which(res[[1]] == 1)
# variables chosen on second round
which(res[[1]] == 2)
# plug in error rates for each round
res[[2]]

# plot y against fitted value from round 1
plot(y, res[[3]][, 1])
# plot y against fitted value from round 2
plot(y, res[[3]][, 2])

# add cross validation
# 2 rounds 10 fold cross validation
res <- DeleteRepeat(x, y, method = HGmultc, xvalid = TRUE, rounds = 2,
    fold = 10, trace = TRUE)

# variables chosen on first round
which(res[[1]] == 1)
# variables chosen on second round
which(res[[1]] == 2)
# cross validated error rates for each round
res[[2]]

# plot y against cross validated fitted value from round 1
plot(y, res[[3]][, 1])
# plot y against cross validated fitted value from round 2
plot(y, res[[3]][, 2])
}

\keyword{ models }
\keyword{multivariate}
