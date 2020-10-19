\name{completeGene}
\alias{completeGene}

\title{Delete genes with any NAs}

\description{
  Removes columns (genes) that contain NAs
}

\usage{
completeGene(Xmat)
}

\arguments{
  \item{Xmat}{A data matrix of cDNA or Affymetrix chip data. The rows are arrays
    and the columns genes.}
}

\value{
  Data matrix with columns containing NAs removed.
}

\author{Glenn Stone}

\seealso{\code{\link{completeArray}}, \code{\link{impute}}}

\keyword{manip}
