\name{HGengine}
\alias{HGengine}

\title{User defined GeneRave models}
\description{
Fits a user defined  model using sparsity priors
}
\usage{
HGengine(X, y, event, weights = rep(1, nrow(X)), bbess = 1e+07,
    kbess = 0, scale, bhat, Le, DLde, D2Lde, HGsc,
    sparsity.prior = "NG", no.prior = NULL,
    Offset = NULL, Offset.par.update = NULL, a = NULL,
    control = HGcontrol())
}

\arguments{
	\item{X}{ real matrix with n rows (samples) and p columns (variables or genes) }
	\item{y}{ real n by 1 vector of response values }
	\item{event}{a vector required for survival models, however could be used as a working vector for other models.]
              Usually set to NULL. }
	\item{weights}{real n by 1 vector of observation weights, default to 1 }
	\item{bbess}{ scale parameter in prior - range 1 to 1e7, default 1e7 }
	\item{kbess}{ shape parameter in prior - 0 to 1, default 0 }
	\item{scale}{scale parameter(s) in the model. If negative, scale
		fixed to abs(scale), otherwise scale estimated using scale as initial
		value.}
	\item{bhat}{initial estimate of beta vector }
	\item{Le}{ function to compute the likelihood function}
	\item{DLde}{function to compute the first derivative of the likelihood function }
	\item{D2Lde}{function to compute the second derivative of the likelihood function  }
	\item{HGsc}{function to update the scale parameter(s) in the model}
	\item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R" }
	\item{no.prior}{ Which terms should not be penalised }
	\item{Offset}{function to compute offset to be added to linear predictor}
	\item{Offset.par.update}{function to update any parameters in the offset }
	\item{a}{intial values of the parameters in the offset}
	\item{control}{ See \code{\link{HGcontrol}} }
}

\details{Fits the user defined model using the EM algorithm described in Kiiveri(2003).
The user must define the following functions
\describe{
\item{Le}{ a function of eta, y, event, weights, scale to compute the likelihood
          function as a function of eta}
\item{DLde}{ a function of eta, y, event, weights, scale to compute the
            first derivative of the likelihood as function of eta }
\item{D2Lde}{ a function of eta, y, event, weights, scale to compute the
             second derivative of the likelihood as fn of eta, or the expected second
             derivative, or the diagonal of some negative (semi) definite matrix.
	           Can be a n by 1 vector (diags) or an n by n neg def matrix}
\item{HGsc}{ a function of eta, y, event, weights, scale for updating the scale parameter(s)}
\item{Offset}{ a function of event, y, a, to compute the offset to add to the linear predictor}
\item{Offset.par.update}{ a function of event, y, a, eta, weights, to update the parameters
                        in the offset function}
}

The sparsity priors which have been implemented are:

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
An object of class HGengine with components:
\item{beta}{p by 1 vector of parameter estimates}
\item{S}{p by 1 logical vector with value TRUE if variable selected i.e
 the corresponding beta is nonzero }
\item{eta}{the estimated linear predictor in the model, eta = X \%*\% beta}
\item{scale}{estimated value of scale parameter(s)}
\item{a}{estimated values of parameter(s) in the offset}
}

\references{
Fan, J. and Li, R. (2001) Variable selection via penalized likelihood.
Journal of American Statistical Association, 96,1348-1360.

Figueiredo, M. (2003) Adaptive sparseness for supervised learning, IEEE
Transactions on Pattern Analysis and Machine Intelligence - PAMI, vol. 25,
no. 9 pp. 1150-1159

Kiiveri, H.T. (2003) A Bayesian approach to variable selection when the number of
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
\code{\link{HGcontrol}}
}

\examples{
# Implement linear regression model using HGengine
# generate data
x <- matrix(rnorm(200 * 200), ncol = 200)
# add intercept to x
X<-cbind(rep(1,200),x)
y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)
# define  functions required by HGengine
Le.g <- function(eta,y,event,weights,scale)
{
	n <- sum(weights)
	L <- -n*log(scale)/2 -0.5*sum(weights*(y-eta)^2)/scale
	L
}

DLde.g <- function(eta,y,event,weights,scale)
{
	weights*(y-eta)/scale
}

D2Lde.g <- function(eta,y,event,weights,scale)
{
	-weights/scale
}

HGsc.g <- function(eta,y,event,weights,scale)
{
	scale <- sum(weights*(y-eta)^2)/sum(weights)
	scale <- max(scale,1e-4)
	scale
}
# generate initial value
bhat<-rnorm(201)  # not optimal but OK for example
# fit sparse linear regression model with scale fixed to 1 and intercept term not penalised
res<-HGengine(X,y,event=NULL, scale=-1,bhat=bhat,Le=Le.g,DLde=DLde.g,D2Lde=D2Lde.g,
              HGsc=HGsc.g,no.prior=1)
#look at estimated coefficients
res$beta[res$S]
# which variables selected
which(res$S)-1
# plot y against fitted value
plot(res$eta , y)
}

\keyword{ misc }
