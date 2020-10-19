\name{eval.surv.fit}
\alias{eval.surv.fit}
\title{Evaluate goodness of fit of a survival model}
\description{
 Compute expected and observed numbers of events in a specified set of intervals given
 a set of survival curves
}
\usage{
eval.surv.fit(times, event, surv, breaks = quantile(times, probs = c(.25, .5, .75)))
}
\arguments{
  \item{times}{real n by 1 vector of survival times  }
  \item{event}{ real n by 1 censoring indicator. The value one denoting corresponding
                survival time uncensored, 0 censored.  }
  \item{surv}{an n by n matrix of survival curves with rows corresponding to times }
  \item{breaks}{ a set of break points to divide the interval from 0 to max(times)
                 but not including 0 or max(times) }
}
\details{
 The range of survival times is divided into intervals by using the breaks vector.
 The observed number of events occuring in each interval is calculated. Censored observations contribute
 the probability of the event occuring in any of the specified intervals using the individuals survival curve.
 The expected number of events occurring in any interval is calculated from each individuals survival curve
 by computing the probability of an event occuring in each interval. This is experimental use with caution.
}
\value{
  \item{observed }{The observed number of events in the specified intervals}
  \item{expected }{The expected number of events in the specified intervals}
  \item{cut.times}{The vector of break points used to define the intervals}
}
\author{ Harri Kiiveri}
\note{ Cross validated survival curves can be used in this function
}
\seealso{ HGcoxreg, HGsurv, xval.survival }
\examples{
set.seed(123)
# generate data
x <- matrix(rnorm(50 * 200), nrow = 50)
lp <- 3 + x[, 1]

# weibull shape parameter
a <- 10

# simulate survival times
time <- rweibull(50, a, exp(-lp / a))

# no censoring
event <- rep(1, 50)

# get cross validated survival curves
res <- xval.survival(x, time, event, fold = 10, trace = TRUE)

# compute observed and fitted counts in default intervals
tmp <- eval.surv.fit(time, event, res$xval.surv)

# observed counts
tmp$observed

# expected counts
tmp$expected
}

\keyword{ survival }
