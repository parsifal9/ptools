\name{HGmultc}
\alias{HGmultc}

\title{GeneRave multiclass classification }

\description{
  Fits multiclass logistic regression models with sparsity priors
}

\usage{
HGmultc(X, y, weights = rep(1, nrow(X)),initb=NULL, sparsity.prior = "NG", lambda = 1,
	bbess = 1e+07, kbess = 0, adj.prior = TRUE, control = HGcontrol())
}

\arguments{
  \item{X}{real matrix with n rows (samples) and p columns (variables or genes)}
  \item{y}{real n by 1 vector of class labels (1, 2, \ldots, C)}
  \item{weights}{real n by 1 vector of observation weights }
  \item{initb}{real (p+1) by G matrix of intial values of class parameters}
  \item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R"}
  \item{lambda}{ridge regression parameter used in computing initial value. Should
   only need modification if \code{X} is illconditioned }
  \item{bbess}{ scale parameter in prior - range 1 to 1e7- default 1e7 }
  \item{kbess}{ shape parameter in prior - range 0 to 1 - default 0}
  \item{adj.prior}{ flag for prior adjustment}
  \item{control}{ See \code{\link{HGcontrol}} }
}

\details{
A symmetric (redundant) parameterisation is used and the prior adjusted
accordingly if \code{adj.prior=TRUE}. If \code{adj.prior=FALSE}, no prior
adjustment is made. No adjustment may be useful when the number of classes is
small, eg 2.

The intercept is forced into the model for each class and is not subject to
the prior.

The default values for \code{bbess} and \code{kbess} give Jeffreys hyperprior.

The sparsity priors implemented are
\describe{
	\item{NG}{ the Normal Gamma prior with gamma scale and shape parameters
		\code{bbess} and \code{kbess} respectively. This default prior appears
		to work well in practise.
	}
	\item{NEG}{ the normal exponential gamma prior with gamma scale and shape
		parameters \code{bbess^-2} and \code{kbess} respectively.
		Griffin and Brown (2005)
	}
	\item{DEG}{ the double exponential gamma prior with gamma scale and shape
		parameters \code{bbess} and \code{kbess} respectively,
		Cawley and Talbot (2006). A special case of the NG prior
	}
	\item{SCAD}{ the SCAD prior (Fan and Li (2001)). The \code{a} and
		\code{lambda} parameters defined by Fan and Li are \code{bbess} and
		\code{kbess} respectively. The values \code{bbess = 3.7} and
		\code{kbess = 2} are Fan and Li's defaults. Experimental use with caution.
	}
	\item{R}{the \dQuote{ridge} prior - proportional to exp(-bbess*beta'beta),
		\code{kbess} is ignored. Does not produce sparse models but may
		sometimes be useful for initial value generation
	}
}

	Note: the default values for \code{bbess} and \code{kbess} assume the
	\dQuote{NG} prior. These values should be specifically set for the other
	priors.
}

\value{
An object of class HGmultc, which is a list with the following components:
\item{beta}{ p+1 by C matrix of parameter estimates - one column for each class}
\item{S}{p+1 by C logical matrix with value true if variable selected i.e
	the corresponding beta is non zero}
\item{P}{n by C matrix of fitted probabilities}
\item{class}{n by 1 vector of fitted class labels}
}

\references{
Fan, J. and Li, R. (2001) Variable selection via penalized likelihood.
Journal of American Statistical Association, 96,1348-1360.

Figueiredo, M. (2003) Adaptive sparseness for supervised learning, IEEE
Transactions on Pattern Analysis and Machine Intelligence - PAMI, vol. 25,
no. 9 pp. 1150-1159

Kiiveri, H.T., A Bayesian approach to variable selection when the number of
variables is very large. In: Science and Statistics: A festschrift for Terry
Speed. IMS Lecture Notes - Monograph Series, Volume TBD (40 or 41).

Kiiveri, H.T. (2008). A general approach to simultaneous model fitting and variable elimination
in response models for biological data with many more variables than observations.
BMC Bioinformatics 2008, 9:195 (15 Apr 2008). Available at
\url{http://www.biomedcentral.com/1471-2105/9/195}

Griffin, J. E., and Brown, P. J. (2005). Alternative prior distributions
for variable selection with very many more variables than observations.
Technical report available at
\url{http://www2.warwick.ac.uk/fac/sci/statistics/crism/research/2005/paper05-10/05-10w.pdf}

Cawley, G., and Talbot, N.,(2006)
Gene selection in cancer classification using sparse logistic
regression with Bayesian regularization
Bioinformatics 22,19,2348-2355.
}

\section{Warning}{
A column of ones is added to the beginning of the \code{X} matrix internally
so variable identifiers are incremented by one. See examples below.
}

\note{
The maximum size matrix for inversion in the algorithm is q by q, where
q = min(n, p). Missing values are not allowed. For models with very large numbers of variables
it may be necesssary to set tolc in HGcontrol to a smaller value such as 1e-4 or 1e-5.
}

\author{ Harri Kiiveri and Glenn Stone }

\seealso{
\code{\link{HGglm}},
\code{\link{HGordcat}},
\code{\link{HGglm}},
\code{\link{HGgaussian}},
\code{\link{HGcoxreg}},
\code{\link{HGsurv}},
\code{\link{DeleteRepeat}}
\code{\link{xvalidate}}
}

\examples{
w <- abs(sin(1:200))
x <- matrix(rnorm(200 * 200), nrow = 200)
y <- (x[, 1] > 0) + 1
k <- 2

# Use default weights and default Jeffreys hyperprior
res <- HGmultc(x, y)

# Display nonzero parameter estimates for class k.
# Note first term is intercept
res$beta[res$S[, k], k]

# Display correct variable identifiers, 0 corresponds to intercept
which(res$S[, k]) - 1

# Display misclassification table
table(y, res$class)

# Use Lasso prior with default scale parameter
res <- HGmultc(x, y, kbess = 1)
table(y, res$class)

# Use weights in w
res <- HGmultc(x, y, weights = w)
table(y,res$class)

# Display matrix of fitted probabilities
round(res$P, 2)

# Compute symmetric parameter values - makes no difference for two class case
bhat <- sweep(res$beta, 1, apply(res$beta, 1, mean))
# Identify chosen variables i.e variables chosen for any class.
bhat[apply(res$S, 1, any), ]
}

\keyword{models}
\keyword{multivariate}
