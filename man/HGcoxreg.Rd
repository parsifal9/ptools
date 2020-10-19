\name{HGcoxreg}

\alias{HGcoxreg}

\title{GeneRave Cox proportional hazards models}

\description{
Fits the Cox proportional Hazards model using sparsity priors. Does simultaneous
variable selection and parameter estimation. Can handle the case of less
observations than variables as well as the case of many more variables than
observations.
}

\usage{
HGcoxreg(x, time, event, weights = rep(1, nrow(x)), sparsity.prior = "NG",
         bbess = 1e7, kbess = 0, b0sc = 1, scale = -1, initb = FALSE,
         no.prior = NULL)
}

\arguments{
	\item{x}{ real matrix with n rows (samples) and p columns (variables or genes)   }
	\item{time}{ real n by 1 vector of survival times  }
	\item{event}{ real n by 1 censoring indicator. One denoting corresponding survival time uncensored, 0 censored.  }
	\item{weights}{ Currently ignored. Real n by 1 vector of observation weights, default to 1 }
	\item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R"}
	\item{bbess}{real number in the range [0, 1e7]}
	\item{kbess}{real number in the range [0, 1]. The value for the SCAD prior can
  be bigger than one, see below.}
	\item{b0sc}{ Length of initial beta vector, defaults to 1. If less than zero use
  internally generated intial value without scaling its length }
	\item{scale}{ Currently ignored }
	\item{initb}{if FALSE generate inital value automatically, if set equal to a p+1 by 1
		vector use this vector to initialise iterations. If \code{initb} is set
		\code{b0sc} is ignored.}
	\item{no.prior}{ Default is to penalise everything }
}

\details{
Fits the Cox proportional hazards model using \code{x} as the matrix of
potential covariates. Ties are handled by adding a small perturbation to
equal survival times. The intercept is NOT forced into the model.

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
		\code{kbess = 2} are Fan and Li's defaults.
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
An object of class HGcoxreg which is a list with components:
\item{beta}{p by 1 vector of parameter estimates}
\item{S}{
	p by 1 logical vector with value true if variable selected
	i.e the corresponding beta is nonzero
}
\item{lp}{
	an n by 1 vector contaning the fitted linear predictor for the model
}
\item{surv0}{
	an n by 3 matrix with columns time, baseline survival function, event
}
\item{survvall}{
	an n by (n + 1) matrix with columns t, estimated survival function for
	subject 1, estimated survival function for subject 2, \ldots, estimated
	survival function for subject n
}
\item{varids}{identifiers of variables selected, 0 denotes intercept}
}

\references{
Cox, D.R., and Oakes, D. (1984), Analysis of Survival Data. Chapman & Hall,
London.

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

\author{ Albert Trajstman and Harri Kiiveri }

\seealso{
\code{\link{HGcoxreg}},
\code{\link{HGgaussian}},
\code{\link{HGglm}},
\code{\link{HGmultc}},
\code{\link{HGordcat}},
\code{\link{HGsurv}},
\code{\link{DeleteRepeat}},
\code{\link{xvalidate}}
}

\note{
Missing values are not allowed. The maximum size matrix for inversion
in the algorithm is q by q, where q = min(n, p).
}

\section{Warning}{There is no type checking of parameters to the function}

\examples{
set.seed(123)
# generate data
x <- matrix(rnorm(50 * 200), nrow = 50)
lp <- -3 + x[, 1]

# weibull shape parameter 1 - exponential distribution
a <- 1

# simulate survival times
time <- rweibull(50, a, exp(-lp / a))

# no censoring
event <- rep(1, 50)

# fit model
res <- HGcoxreg(x, time, event, b0sc=0.5)

# variables selected 0 denotes intercept
res$varids

# parameter estimates
res$beta[res$S]
}

\keyword{ models }
\keyword{ survival }
