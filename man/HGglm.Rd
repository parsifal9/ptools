\name{HGglm}

\alias{HGglm}

\title{GeneRave generalised linear models }

\description{
Fits generalised linear models with sparsity priors. Does simultaneous
variable selection and parameter estimation. Can handle the case of fewer
variables than observations as well as the case of many more variables
than observations.
}

\usage{
HGglm(x, f, event = NULL, weights = rep(1, nrow(x)),
      sparsity.prior = "NG", bbess = 1e+07, kbess = 0,
      b0sc = 15, scale = -1, initb = "FALSE", model = "N",
      fvalfn = NULL, ifvalfn = NULL, varfn = NULL, drfn = NULL,
      devfn = NULL, scale.updatefn = NULL, no.prior = 1) }

\arguments{
	\item{x}{real matrix with \eqn{n} rows (samples) and \eqn{p}
	  columns (variables or genes).}
	\item{f}{ real \eqn{n \times 1} vector of response values.}
	\item{event}{ NULL except for model "B" in which case it is the
	  vector of the number of binomial trials for each observation.}
	\item{weights}{real \eqn{n \times 1} vector of observation weights,
	  defaults to 1.}
	\item{sparsity.prior}{one of \code{"NG"}, \code{"NEG"},
	  \code{"DEG"}, \code{"SCAD"}, or \code{"R"}.}
	\item{bbess}{real number in the range \code{[0, 1e+07]}}
	\item{kbess}{real number in the range \code{[0, 1]}. The value
	  for the \code{SCAD} prior can be bigger than one, see below.}
	\item{b0sc}{length of initial beta vector---likely useful
	  values will be in the range 1--25. Ignored for model "B". If
	  less than zero, use internally generated intial value without
	  scaling its length.}
	\item{scale}{value of scale parameter in generalised linear
	  model. If negative, \code{scale} is fixed to
	  \code{abs(scale)}, otherwise scale estimated using
	  \code{scale} as initial value.}
	\item{initb}{if \code{FALSE} generate inital value
	  automatically; if set equal to a \eqn{p + 1 \times 1} vector
	use this vector to initialise iterations. If \code{initb} is
	set, \code{b0sc} is ignored.}
	\item{model}{The possible values for model are:
	  \itemize{
	    \item{\code{N}}{normal disribution}
	    \item{\code{P}}{poisson distribution}
	    \item{\code{B}}{binomial - logistic regression}
	    \item{\code{G}}{gamma}
	    \item{\code{IG}}{inverse gaussian distribution}
	    \item{\code{Own}}{define your own model, in which case the
	      arguments \code{fvalfn}, \code{ifvalfn}, \code{varfn},
	      \code{drfn}, \code{devfn}, and \code{scale.updatefn} must
	      be specified (see below).}
	  }
	}
	\item{fvalfn}{a function of the linear predictor
	  \eqn{\eta} which returns \eqn{\mu}, the fitted value.
	}
	\item{ifvalfn}{a function of the fitted value \eqn{\mu} which
	  returns \eqn{\eta}, the linear predictor.  Hence,
	  \code{ifvalfn} is the inverse function of \code{fvalfn}.}
	\item{varfn}{a function of \eqn{\mu} which returns the value of
	  the variance function}
	\item{drfn}{function of \eqn{\mu} which returns d(mu)/d(eta)}
	\item{devfn}{
	  a function of \code{mu}, \code{y}, \code{event}, \code{scale},
	  and \code{weights} which returns the value of the deviance.
	}
	\item{scale.updatefn}{
	  a function of \code{mu}, \code{y}, \code{event}, \code{scale},
	  and \code{weights} which returns the value of the scale
	  parameter.
	}
	\item{no.prior}{Vector of column numbers specifying which
	  variables are not to be subject to the prior. Default is 1,
	  i.e., the intercept is not penalised. Can be set to
	  \code{NULL} if intercept is to be penalised.}
}

\details{
  The functions supplied to the arguments \code{fvalfn},
  \code{ifvalfn}, \code{varfn}, \code{drfn}, \code{devfn}, and
  \code{scale.updatefn} are only required when \code{model = "Own"}.
  Some experimentation with the inital value generation may be required
  for \code{"Own"} models.

  The intercept is forced into the model by default and is not subject to
  the prior. It is always added to the begining of \code{x}.

  The sparsity priors implemented are:
  \describe{
    \item{NG}{the Normal Gamma prior with gamma scale and shape
      parameters bbess and kbess respectively. This default prior appears
      to work well in practice.}

    \item{NEG}{the normal exponential gamma prior with gamma
      scale and shape parameters \code{bbess}\eqn{^-2} and \code{kbess}
      respectively (Griffin and Brown (2005)).}

    \item{DEG}{the double exponential gamma prior with gamma
      scale and shape parameters \code{bbess} and \code{kbess}
      respectively (Cawley and Talbot (2006)). A special case of the NG
      prior.}

    \item{SCAD}{the SCAD prior (Fan and Li(2001)). The parameters
      \eqn{a} and \eqn{\lambda} defined by Fan and Li are \code{bbess}
      and \code{kbess} respectively. The values \code{bbess = 3.7} and
      \code{kbess = 2} are Fan and Li's defaults.}

    \item{R}{the ridge prior---proportional to \code{exp(-bbess *sum(beta * beta))},
     where \code{beta} is the vector of parameters. The parameter \code{kbess}
     is ignored. Does not produce sparse models but may sometimes be  useful for initial
      value generation.}
  }
  Note: the default values for \code{bbess} and \code{kbess} are
  sensible values when the \code{NG} prior has been chosen.  These
  values should be specifically set for the other priors.
}

\value{An object of class \code{"HGglm"} that is a list with the
  following components:

 \item{beta}{\eqn{p + 1 \times 1} vector of parameter estimates.}

 \item{S}{\eqn{p + 1 \times 1} logical vector with value true if
   variable selected, i.e., the corresponding beta is nonzero.}

 \item{fv}{\eqn{n \times 1} vector of fitted values.}

 \item{varids}{identifiers of variables selected; 0 denotes intercept.}

 \item{sc}{estimated value of scale parameter.}

 \item{model}{The model used, i.e., one of \code{"N"}, \code{"P"},
   \code{"B"}, \code{"G"}, \code{"IG"}, or \code{"Own"}.}

 \item{bd}{Only for model \code{"B"}. The number of trials in the
   binomial for each sample.}
}

\note{The maximum size matrix for inversion in the algorithm is \eqn{q
    \times q}, where q = min(n, p). Missing values are not
  allowed.
}

\references{
  McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear
    Models}, 2nd ed., Chapman and Hall: London.

  Fan, J. and Li, R. (2001) Variable selection via penalized likelihood.
  \emph{Journal of American Statistical Association}, \textbf{96},
  1348--1360.

  Figueiredo, M. (2003) Adaptive sparseness for supervised
  learning. \emph{IEEE Transactions on Pattern Analysis and Machine
    Intelligence}, \textbf{25} (9), 1150--1159.
    
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

  Cawley, G., and Talbot, N. (2006) Gene selection in cancer
  classification using sparse logistic regression with Bayesian
  regularization. \emph{Bioinformatics}, \textbf{22} (19), 2348--2355.
}

\author{ Harri Kiiveri }

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

\examples{
#
# EXAMPLE 1 -logistic regression using HGglm
#
x <- matrix(rnorm(200 * 200), nrow = 200)
y <- as.numeric(x[, 1] > 0)

# Fit glm with default weights and default Jeffreys hyperprior
res <- HGglm(x, y, model = "B")

# Display nonzero parameter estimates. Note first term is intercept
res$beta[res$S]

# Display identifiers of selected variables - zero corresponds to intercept
res$varids

# Display misclassification table for logistic regression
table(y, res$fv > 0.5)

# Use Lasso prior with default scale parameters
res <- HGglm(x, y, model = "B", kbess = 1)
res$beta[res$S]

# Use weights in w
w <- abs(sin(1:200))
res <- HGglm(x, y, model = "B", weights = w)
res$beta[res$S]

# Plot fitted values against y
plot(y, res$fv)

#
# EXAMPLE 2 - GLM - linear regression using "Own" function
#

# Define required functions for linear regression

fvalfn <- function(eta) { eta }

ifvalfn <- function(mu) { mu }

varfn <- function(mu) { 1 }

drfn <- function(mu) { 1 }

devfn <- function(mu, y, event, scale, weights)
{
	L <- -0.5 * sum(weights) * log(scale) - 0.5 *
		sum(((weights * (y - mu)^2) / scale))
	L
}

scupdate<-function(mu, y, event, scale, weights)
{
	scale <- sum(weights * (y - mu)^2) / length(y)
	scale <- max(scale, 1e-2)
	scale
}

#
# generate data and test
x <- matrix(rnorm(200 * 200), ncol = 200)
y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)

# scale fixed to 1
\dontrun{
res <- HGglm(x, y, model = "Own",
	fvalfn = fvalfn, varfn = varfn, drfn = drfn, devfn = devfn,
	scale.updatefn = scupdate)
}

# parameter estimates
res$beta[res$S]
res$varids

# also estimate scale parameter
\dontrun{
res <- HGglm(x, y, scale = 1, b0sc = 1, model = "Own",
	fvalfn = fvalfn, ifvalfn = ifvalfn, varfn = varfn, drfn = drfn, devfn = devfn,
	scale.updatefn = scupdate)
}

# parameter estimates
res$beta[res$S]

# variable ids 0 is intercept
res$varids

# estimated standard deviation
sqrt(res$scale)

# restart iterations from previous values
\dontrun{
res <- HGglm(x, y, scale = 1, b0sc = 1, initb = res$beta, model = "Own",
	fvalfn = fvalfn, ifvalfn = ifvalfn, varfn = varfn, drfn = drfn, devfn = devfn,
	scale.updatefn = scupdate)
}
#
# EXAMPLE 3- GLM - logistic regression using "Own" function
#
x <- matrix(rnorm(200 * 200), nrow = 200)
y <- as.numeric(x[, 1] > 0)
#
cupcap <- function (x, lim)
{
	pmin(pmax(x, -lim), lim)
}

fvalfn <- function (eta)
{
	eta <- cupcap(eta, 100)
	exp(eta) / (1 + exp(eta))
}

devfn <- function (mu, y, event, scale, weights)
{
	mu <- pmin(1 - 1e-10, mu)
	mu <- pmax(1e-10, mu)
	sum(weights * y * log(mu) + weights * (1 - y) * log(1 - mu))
}

drfn <- function (mu)
{
	mu <- pmin(1 - 1e-10, mu)
	mu <- pmax(1e-10, mu)
	mu * (1 - mu)
}

varfn <- function (mu)
{
	mu <- pmin(1 - 1e-10, mu)
	mu <- pmax(1e-10, mu)
	mu * (1 - mu)
}

scupdate <- function (mu, y, event, scale, weights)
{
	scale
}

res <- HGglm(x, y, scale = 1, b0sc = 1, initb = rnorm(201), model = "Own",
	fvalfn = fvalfn, varfn = varfn, drfn = drfn, devfn = devfn,
	scale.updatefn = scupdate)

} % end examples

\keyword{ models  }
\keyword{ regression }
