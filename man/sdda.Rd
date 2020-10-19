\name{sdda}
\alias{sdda}
\alias{sdda.matrix}
\alias{sdda.dda}

\title{Stepwise Diagonal Discriminant Analysis}
\description{
  Uses a forwards stepwise strategy to add variables one at a time to a
  diagonal discriminant rule.
}
\usage{
sdda(obj, \dots)
sdda.matrix(obj, y, priors, \dots)
sdda.dda(obj, X, y, start = rep(FALSE, p), never = rep(FALSE, p),
    useprob = TRUE, usecache = TRUE, usexval = TRUE, jmax= -20,
    \dots)
}

\arguments{
	\item{obj}{An object of class \code{\link{dlda}} or \code{\link{dqda}}}
	\item{X}{Training data matrix - rows are observations columns of variables.}
	\item{y}{Vector of true class labels}
	\item{priors}{ ~ TODO ~ }
	\item{start}{
		Logical vector indicating which variables (if any) should start in the
		rule.
	}
	\item{never}{
		Logical vector indicating which variables (if any) should never be in
		the rule.
	}
	\item{useprob}{ ~ TODO ~ }
	\item{usecache}{ ~ TODO ~ }
	\item{usexval}{ ~ TODO ~ }
	\item{jmax}{ ~ TODO ~ }
	\item{\dots}{ ~ TODO ~ }
}
\details{
  Uses a forwards stepwise strategy to add variables one at a time to a
  diagonal discriminant rule. At each step the variable that most
  decreases the (leave-one-out) cross-validated error rate is
  added. Adding stops when the cross-validated error rate cannot be decreased.
}
\value{
  \item{}{A copy of \code{obj} plus \ldots}
  \item{S}{Vector of flags of the variables in the rule, including \code{start}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Glenn Stone}

\seealso{\code{\link{dlda}}, \code{\link{dqda}},
  \code{\link{which.genes}}, \code{\link{predict.sdda}}}

%\examples{}
\keyword{models}
\keyword{multivariate}
