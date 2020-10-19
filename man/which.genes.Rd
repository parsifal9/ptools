\name{which.genes}
\alias{which.genes}
\alias{which.genes.fastRP}
\alias{which.genes.HGcoxreg}
\alias{which.genes.HGgaussian}
\alias{which.genes.HGglm}
\alias{which.genes.HGmultc}
\alias{which.genes.HGordcat}
\alias{which.genes.HGsurv}
\alias{which.genes.sdda}

\title{which genes}

\description{
Returns a list of genes (variables) used in the model. Takes care of
intercepts. Methods exist for \code{\link{fastRP}}, \code{\link{HGcoxreg}},
\code{\link{HGgaussian}}, \code{\link{HGglm}}, \code{\link{HGmultc}},
\code{\link{HGordcat}}, \code{\link{HGsurv}} and \code{\link{sdda}}.
}

\usage{
which.genes(object, \dots)
}

\arguments{
\item{object}{
	An object of the appropriate class such as \code{\link{HGmultc}}
}
\item{\dots}{Any other parameters}
}

\value{
}

\seealso{
\code{\link{fastRP}},
\code{\link{HGcoxreg}},
\code{\link{HGgaussian}},
\code{\link{HGglm}},
\code{\link{HGmultc}},
\code{\link{HGordcat}},
\code{\link{HGsurv}},
\code{\link{sdda}}
}

\examples{
x <- matrix(rnorm(200 * 200), nrow = 200)
y <- (x[, 10] > 0) + 1

# Use default weights and default Jeffreys hyperprior
res <- HGmultc(x, y)

which.genes(res)
}

\author{Harri Kiiveri}
\keyword{misc}
