\name{dqda}
\alias{dqda}

\title{Diagonal Quadratic Discriminant Analysis}
\description{
Quadratic Discriminant Analysis with a diagonal covariance matrix.
}

\usage{
dqda(X, y, priors)
}

\arguments{
  \item{X}{Training data matrix - rows are observations columns of variables.}
  \item{y}{Vector of true class labels}
  \item{priors}{Vector of class priors (defaults to equal)}
}

\value{
  \item{means}{matrix of means of variables for each class}
  \item{vars}{matrix of variances of variables for each class}
  \item{counts}{vector of counts of obs. in each class}
}

\author{Glenn Stone}
\note{
Although it can be used on its own, this function is primarily a
service function for use with \code{\link{sdda}}.
}

\seealso{\code{\link{sdda}}, \code{\link{dlda}}, \code{\link{predict.dqda}}}

\keyword{models}
\keyword{multivariate}
