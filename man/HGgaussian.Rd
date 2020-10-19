\name{HGgaussian}

\alias{HGgaussian}

\title{GeneRave linear regression  }

\description{
Fits linear regression models with sparsity priors. Does simultaneous
variable selection and parameter estimation. Can handle the case of less observations than
variables as well as the case of  many more variables than observations.
}

\usage{
HGgaussian(x, f, weights = rep(1, nrow(x)), sparsity.prior="NG",
           bbess=1e7, kbess=0, b0sc = 5, scale = -1, initb = NULL, no.prior = 1,
           control = HGcontrol(eta.lim = 1e30))
}

\arguments{
	\item{x}{ real matrix with n rows (samples) and p columns (variables or genes)  }
	\item{f}{ real n by 1 vector of response values }
	\item{weights}{real n by 1 vector of observation weights, default to 1 }
	\item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R"}
	\item{bbess}{real number in the range [0,1e7]. When kbess equals one for the NG prior
  the usual penalty parameter in L1 regression is (2/bbess)^0.5. In this case values
  of bbess in the interval (0,1] will often be useful}
	\item{kbess}{real number in the range [0,1]. The value for the SCAD prior can
  be bigger than one, see below.}
	\item{b0sc}{scale factor - scales the model to make the largest component of the automatically generated initial beta
    vector equal to b0sc in absolute value - likely useful values in the range 1 - 25. If
    less than zero use internally generated intial value without scaling its length.
		This parameter is ignored if the number of rows of x exceeds the number of columns in which case
    the automatically generated initial value is approximately equal to the usual regression estimate.}
	\item{scale}{value of scale (variance) parameter in the linear regression model. If negative, scale
		fixed to abs(scale), otherwise scale estimated using scale as initial
		value.}
	\item{initb}{if NULL generate inital value automatically, if set equal to a p+1 by 1
		vector use this vector to initialise iterations. If initb is set b0sc is ignored.}
	\item{no.prior}{ Default is 1, don't penalise the intercept.}
	\item{control}{See HGcontrol}
}

\details{
The error distribution is assumed to be Gaussian. The intercept is forced into
the model by default and is not subject to the prior. It is added to the begining of
\code{x}.

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
An object of class HGgaussian which is a list with the following components:
\item{beta}{ p+1 by 1 vector of parameter estimates   }
\item{S}{
	p+1 by 1 logical vector with value true if variable selected i.e
	the corresponding beta is nonzero
}
\item{fv}{n by 1 vector of fitted values  }
\item{varids}{identifiers of variables selected, 0 denotes intercept }
\item{sigma2}{estimated value of error variance parameter  }
}

\note{
The maximum size matrix for inversion in the algorithm is q by q, where
q = min(n, p). Missing values are not allowed. Strange behaviour may result if the
x matrix is ill-conditioned. Try the matrix scale(x) instead.}

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

\author{ Harri Kiiveri }

\seealso{
\code{\link{HGcoxreg}},
\code{\link{HGgaussian}},
\code{\link{HGglm}},
\code{\link{HGmultc}},
\code{\link{HGordcat}},
\code{\link{HGsurv}},
\code{\link{DeleteRepeat}}
\code{\link{xvalidate}}
}

\examples{
# EXAMPLE 1 - generate regression data and test

x <- matrix(rnorm(200 * 200), ncol = 200)
y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)

# Scale fixed to 1
res <- HGgaussian(x, y)

# Parameter estimates
res$beta[res$S]
res$varids

# Slso estimate scale parameter
res <- HGgaussian(x, y, sc = 1, b0sc = 1)

# Parameter estimates
res$beta[res$S]
# Variable ids 0 is intercept
res$varids
# Estimated standard deviation
sqrt(res$sigma2)

# Restart iterations from previous values
res <- HGgaussian(x, y, sc = 1, b0sc = 1, initb = res$beta)

#EXAMPLE 2 - regression spline example

x <- seq(-2, 2, .02) + 1e-8
y <- sin(3 * x) / (3 * x)
# Compute test function and data
yy <- y + 0.125 * rnorm(201)

tmp <- outer(x, x, FUN = "-")
tmp <- tmp^3 * as.numeric(tmp > 0)
# Compute regression spline basis
xx <- cbind(rep(1, 201), x, x^2, x^3, tmp)

res <- HGgaussian(xx, yy)

# Identify chosen knot points
x[res$varids]

# Plot results
plot(x, res$fv)
lines(x, y, col = "red")
points(x, yy, col = "green")
}

\keyword{ models }
\keyword{ regression }
