#Reference
#ftp://whitechapel.media.mit.edu/pub/tech-reports/TR-514-ABSTRACT.html
#select no of factors for profile analysis

bayesk <- function(l, N, d, k)
{
### l is eigenvalues of matrix = singular values squared
### N is no of variables(genes)
### k is number of basis functions being considered
### d is rank of the residual matrix
  m <- d * k - k * (k + 1) / 2
  i1 <- 1:k
  i2 <- (k + 1):d
  d1 <- (d - i1 + 1) / 2

  logpu <- -k * log(2) + sum(log(gamma(d1)) - d1 * log(pi))
  v <- sum(l[i2]) / (d - k)# k must be less than d

  lhat <- rep(0, d)
  lhat[i1] <- l[i1]
  lhat[i2] <- v

  logpdk <- logpu - 0.5 * N * sum(log(l[i1])) - 0.5 * N * (d - k) * log(v)
  logpdk <- logpdk + 0.5 * (m + k) * log(2 * pi)

  logdetaz <- 0
  for (i in 1:k) {
    j1 <- (i + 1):d
    logdetaz <- logdetaz + sum(log(1 / lhat[j1] - 1 / lhat[i]) +
    	log(l[i] - l[j1]) + log(N))
  }
  logpdk <- logpdk - 0.5 * logdetaz - 0.5 * k * log(N)

  logpdk
}
