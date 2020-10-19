\name{impute}
\alias{impute}

\title{Simple imputation of Microarray data}

\description{
  Imputes missing values in microarray data matrix using a simple row
  mean + column mean - grand mean method.
}

\usage{
impute(Xmat)
}

\arguments{
  \item{Xmat}{A data matrix of cDNA or Affymetrix chip data. The rows are arrays
    and the columns genes.}
}

\value{
  Imputed data matrix.
}

\author{Glenn Stone}

\seealso{\code{\link{completeArray}}, \code{\link{completeGene}}}

\keyword{manip}
