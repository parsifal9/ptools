\name{HGL1reg}
\alias{HGL1reg}
\title{ Fit L1 regression model with sparsity priors}
\description{
Fits L1 linear regression models with sparsity priors. Does simultaneous
variable selection and parameter estimation. Can handle the case of less observations than
variables as well as the case of many more variables than observations.
}
\usage{
HGL1reg(x, f, weights = rep(1, nrow(x)), sparsity.prior = "NG", bbess = 1e+07, kbess = 1, b0sc = 5, scale = -1, initb = NULL, no.prior = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ real matrix with n rows (samples) and p columns (variables or genes)  }
	\item{f}{ real n by 1 vector of response values }
	\item{weights}{real n by 1 vector of observation weights, default to 1 }
	\item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R"}
	\item{bbess}{real number in the range [0,1e7]. When kbess equals one for the NG prior
  the usual penalty parameter in L1 regression is (2/bbess)^0.5. In this case values
  of bbess in the interval (0,1] will often be useful}
	\item{kbess}{real number in the range [0,1]. The value for the SCAD prior can
  be bigger than one, see for example the help for HGmultc.}
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
}
\details{
The error distribution is assumed to be a double exponential and is formulated as the mixture
of a normal distribution and a gamma distribution. The adjustment to prevent zero divisors given
in the reference below is used. The intercept is forced into
the model by default and is not subject to the prior. It is added to the begining of
\code{x}
}
\value{
An object of class HGL1reg which is a list with the following components:
\item{beta}{ p+1 by 1 vector of parameter estimates   }
\item{S}{
	p+1 by 1 logical vector with value true if variable selected i.e
	the corresponding beta is nonzero
}
\item{fv}{n by 1 vector of fitted values  }
\item{varids}{identifiers of variables selected, 0 denotes intercept }
\item{sigma2}{estimated value of error variance parameter  }
}
\references{Lange, K. L. and Sinsheimer, J. S. (1993), Normal/Independent Distributions
and Their Applications in Robust Regression, Journal of Computational and Graphical Statistics, 2, 175-198 }
\author{Harri Kiiveri }
\note{ The maximum size matrix for inversion in the algorithm is q by q, where
q = min(n, p). Missing values are not allowed.
}
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
# EXAMPLE  - generate regression data and test
x <- matrix(rnorm(200 * 200), ncol = 200)
y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)
# Scale fixed to 1
res <- HGL1reg(x, y,kbess=0,bbess=0.1)
# Parameter estimates
res$beta[res$S]
res$varids
# get fitted values
fv<-predict(res,x)
plot(y,fv)
# do tenfold cross validation - not run
# tmp<-xvalidate(x,y,method=HGL1reg,kbess=0.5,bbess=0.1,fold=10,trace=T)
# plot(y,tmp)
}

\keyword{ models}
\keyword{ multivariate }
