\name{completeArray}
\alias{completeArray}

\title{Delete Microarrays with any NAs}

\description{
  Removes rows (microarrays) that contain NAs
}

\usage{
completeArray(Xmat)
}

\arguments{
  \item{Xmat}{A data matrix of cDNA or Affymetrix chip data. The rows are arrays
    and the columns genes.}
}

\value{
  Data matrix with rows containing NAs removed.
}

\author{Glenn Stone}

\seealso{\code{\link{completeGene}}, \code{\link{impute}}}

\keyword{manip}
