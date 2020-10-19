\name{predict.dlda}
\alias{predict.dlda}

\title{Predict method for diagonal discriminant analysis}
\description{
Produces predicted classes or class probabilities for a \code{\link{dlda}}
rule.
}

\usage{
\method{predict}{dlda}(object, newdata, type = c("class", "prob"), \dots)
}

\arguments{
\item{object}{A fitted \code{dlda} rule object.}
\item{newdata}{
	A dataset to make predictions for. Unlike many predict methods, this is
	required, as \code{\link{sdda}} does not store the original data.
}
\item{type}{
	If \dQuote{class} return the predicted classes, else if \dQuote{prob}
	the class probabilities - note \dQuote{class} is the default.
}
\item{\dots}{
parameters passed to the base \code{\link{dlda}} or \code{\link{dqda}}
predict method. See \code{type} in \code{\link{predict.dlda}} and
\code{\link{predict.dqda}}.
}
}

\value{
If \code{type} is \dQuote{class} a vector of predicted class labels, else if
\code{type} is \dQuote{prob} then a matrix of class probabilities, each row
representing the probabilities for an observation in \code{newdata}.
}

\author{Glenn Stone}

\seealso{\code{\link{predict.sdda}} and \code{\link{predict.dqda}} }

\keyword{models}
\keyword{multivariate}
