\name{break.ties}
\alias{break.ties}
\title{Survival time tie breaker}
\description{Breaks ties in a vector of survival times

}
\usage{
break.ties(x, sd.divisor = 10)
}

\arguments{
  \item{x}{vector of survival times}
  \item{sd.divisor} {standard deviation divisor}
}

\details{
Adds a random normal variate to tied survival times. The normal variate has a standard deviation
equal to the minimum distance between non tied observations divided by sd.divisor.
}

\value{ a vector of survival times with ties broken }

\author{Albert Trajstman}

\seealso{ \code{\link{xval.survival}} }

\examples{
times <- rnorm(10)
times[5] <- times[1]
cat(times, "\n")
break.ties(times)
}

\keyword{misc}
