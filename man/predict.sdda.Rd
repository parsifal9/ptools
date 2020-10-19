\name{predict.sdda}
\alias{predict.sdda}

\title{Predict method for stepwise diagonal discriminant analysis}
\description{
Produces predicted classes or class probabilities for an \code{\link{sdda}}
rule.
}

\usage{
\method{predict}{sdda}(object, newdata, \dots)
}

\arguments{
\item{object}{A fitted \code{\link{sdda}} rule object.}
\item{newdata}{
A dataset to make predictions for. Unlike many predict methods, this is
required, as \code{\link{sdda}} does not store the original data.
}
\item{\dots}{
parameters passed to the based \code{dlda} or \code{dqda}
predict method. See \code{type} in \code{\link{predict.dlda}} and
\code{\link{predict.dqda}}.}
}

\value{
See \code{\link{predict.dlda}} and \code{\link{predict.dqda}}.
}

\author{Glenn Stone}

\seealso{\code{\link{predict.dlda}} and \code{\link{predict.dqda}}  }

\keyword{models}
\keyword{multivariate}
