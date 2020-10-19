\name{HGsurv}

\alias{HGsurv}

\title{GeneRave parametric proportional hazards models}

\description{

Fits parametric proportional hazards models using a sparsity prior. Does simultaneous
variable selection and parameter estimation. Can handle the case of less observations than
variables as well as the case of  many more variables than observations.
}

\usage{
HGsurv(x, time, event, weights = rep(1, nrow(x)), sparsity.prior="NG",
      bbess=1e7, kbess=0, b0sc = 1, scale = -1, initb = "FALSE",
      loginthaz = NULL, aupdate = NULL, a = NULL, no.prior = 1)
}

\arguments{
	\item{x}{ real matrix with n rows (samples) and p columns (variables or genes)   }
	\item{time}{ real n by 1 vector of survival times  }
	\item{event}{ real n by 1 censoring indicator. The value one denoting corresponding survival time uncensored, 0 censored.  }
	\item{weights}{ real n by 1 vector of observation weights, default to 1 }
	\item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R"}
	\item{bbess}{real number in the range [0,1e7]}
	\item{kbess}{real number in the range [0,1]. The value for the SCAD prior can
   be bigger than one, see below.}
	\item{b0sc}{ Length of initial beta vector, defaults to 1. If less than zero use
  internally generated intial value without scaling its length}
	\item{scale}{ Currently ignored }
	\item{initb}{if FALSE generate inital value automatically, if set equal to a p+1 by 1
		vector use this vector to initialise iterations. If initb is set b0sc is ignored.}
	\item{loginthaz}{ function for computing log integrated hazard function }
	\item{aupdate}{function to update estmates of parameters in log integrated hazard
		function}
	\item{a}{ intial values for log integrated hazard function parameter(s) a }
	\item{no.prior}{ default is 1, do not penalise the intercept }
}

\details{
Fits parametric proportional hazards model using \code{x} as the matrix of
potential covariates. The model is converted internally to a Poisson
regression model as per Aitkin and Clayton(1980) below.
The intercept is forced into the model and is not subject to the prior. It is
added to the begining of \code{x}.

The following functions must be specified to define the parametric model.
\describe{
\item{loginthaz}{the log integrated hazard function, a function of event, time, and a}
\item{aupdate}{a function for updating shape parameters in the hazard function, a
      function of event, time, a0, eta, and weights}
\item{a}{a vector of initial values for the parameters in the log integrated hazard function}
}
For details see the example below.

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
An object of class HGsurv which is a list with components:
\item{beta}{p+1 by 1 vector of parameter estimates}
\item{S}{
	p+1 by 1 logical vector with value true if variable selected
	i.e the corresponding beta is nonzero
}
\item{lp}{
an n by 1 vector contaning the fitted linear predictor for the model plus the offset.
}
\item{surv0}{
an n by 3 matrix with columns time, baseline survival function (beta = 0),
event
}
\item{survall}{
an n by (n+1) matrix with columns time, estimated survival function for
subject 1, estimated survival function for subject 2, \ldots,
estimated survival function for subject n
}
\item{varids}{identifiers of variables selected, 0 denotes intercept}
\item{haz.par}{estimated (shape) parameter(s) in the hazard function}
\item{fv}{
the fitted value for the poisson model, with \eqn{i^{th}}{ith} component
\eqn{\lambda_{t_i} \times exp(\beta'x[i, ])}{lambda[t[i]] * exp(beta'x[i, ])},
where \eqn{\lambda}{lambda} is the integrated hazard function
}
}

\references{
Aitkin, M. and Clayton, C. (1980), The Fitting of Exponential, Weibull and
Extreme Value Distributions to Complex Censored Survival Data using GLIM.
Applied Statistics, 29: 156-163.

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

\author{Harri Kiiveri}

\seealso{
\code{\link{HGglm}},
\code{\link{HGordcat}},
\code{\link{HGmultc}},
\code{\link{HGgaussian}},
\code{\link{HGcoxreg}},
\code{\link{HGsurv}},
\code{\link{DeleteRepeat}}
\code{\link{xvalidate}}
}

\note{ Missing values are not allowed.The maximum size matrix for inversion in the algorithm
 is q by q where q = min(n, p).}

\section{Warning }{There is no type checking of parameters to the function}

\examples{
# functions required for weibull
wloginthaz <- function (event, time, a)
{
	a * log(time)
}

waupdate <- function (event, time, a0, eta, weights)
{
# a0 is current value of a
	p0 <- exp(eta)
	a1 <- sum(log(time) * (p0 - event) * weights)
	a1 <- sum(weights * event) / a1
	a0 + 0.1 * (a1 - a0)
}

# generate data
set.seed(123)
x <- matrix(rnorm(50 * 200), nrow = 50)
lp <- -3 + x[, 1]

# weibull shape parameter 0.5
a <- 0.5

# simulate survival times
time <- rweibull(50, a, exp(-lp / a))

# no censoring
event <- rep(1, 50)

# fit model
res <- HGsurv(x, time, event, loginthaz = wloginthaz, aupdate = waupdate, a = 0.25)

# variables selected 0 denotes intercept
res$varids

# parameter estimates
res$beta[res$S]
}

\keyword{ models }
\keyword{ survival }
