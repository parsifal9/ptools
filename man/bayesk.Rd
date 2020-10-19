\name{bayesk}
\alias{bayesk}

\title{ function to select number of principal components }

\description{ function to select number of principal components }

\usage{
bayesk(l, N, d, k)
}

\arguments{
  \item{l}{ eigenvalues of matrix = singular values squared }
  \item{N}{ no of variables(genes) }
  \item{d}{ number of basis functions being considered }
  \item{k}{ rank of the residual matrix }
}

\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

\dontrun{
X <- t(X)
n <- dim(X)[1]
k <- dim(X)[2]
k1 <- dim(TT)[1]
r <- dim(TT)[2]

PB <- TT \%*\% solve(t(TT) \%*\% TT) \%*\% t(TT)
H <- diag(rep(1, k)) - matrix(1, k, k) / k
PW <- H \%*\% (diag(rep(1, k)) - PB)
R <- X \%*\% PW / sqrt(k)
tt <- svd(R)

temp <- rep(0,80)
for (i in 1:80) {
	cat(bayesk(tt$d^2, 19968, 44, i), "\n")
	temp[i] <- bayesk(tt$d^2, 19968, 44, i)
	cat(i, " ", temp[i], "\n")
}
plot(temp)
} % end dontrun
} % end examples

\keyword{ misc }
