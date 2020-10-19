\name{dlda}
\alias{dlda}

\title{Diagonal Linear Discriminant Analysis}
\description{
  Linear Discriminant Analysis with a diagonal covariance matrix.
}

\usage{
dlda(X, y, priors)
}

\arguments{
  \item{X}{Training data matrix - rows are observations columns of variables.}
  \item{y}{Vector of true class labels}
  \item{priors}{Vector of class priors (defaults to equal)}
}

\value{
  \item{means}{matrix of means of variables for each class}
  \item{vars}{vector of variances of variables}
  \item{counts}{vector of counts of obs. in each class}
}

\author{Glenn Stone}
\note{
Although it can be used on its own, this function is primarily a
service function for use with \code{\link{sdda}.}
}

\seealso{\code{\link{sdda}}, \code{\link{dqda}}, \code{\link{predict.dlda}}}

%\examples{}
\keyword{models}
\keyword{multivariate}
