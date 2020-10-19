\name{xval.survival}
\alias{xval.survival}
\title{Cross validated survival curves}
\description{
Produces cross validated survival curves and the cross validated linear predictor
for proportional hazard models
}
\usage{
xval.survival(X, y, event , method = HGcoxreg, fold = NULL, trace = FALSE,
              weights = rep(1, nrow(X)), ...)
}
\arguments{
  \item{X}{real matrix with n rows (samples) and p columns (variables or genes)  }
  \item{y}{ real n by 1 vector of survival times  }
  \item{event}{ real n by 1 censoring indicator. The value one denoting corresponding
           survival time uncensored, 0 censored. }
  \item{method}{One of HGcoxreg or HGsurv }
  \item{fold}{an integer specifying the number of folds to use in the cross validation.
	            Default value NULL gives leave one out cross validation}
  \item{trace}{a logical variable with value TRUE if information about cross validation
	            progress to be output and value FALSE if no output required.}
  \item{weights}{ real n by 1 vector of observation weights, default to 1  }
  \item{\dots}{ any other arguments required by the methods HGcoxreg or HGsurv }
}
\details{
For each group of observations held in, a baseline survival function and coefficient
vector beta is estimated. This baseline survival function is linearly interpolated to give values
for each time in y. The interpolated survival function and coefficient vector beta are then used to
produce survival curves for the observations held out.

}
\value{
A list with components
  \item{xval.lp }{cross validated linear predictor}
  \item{times }{vector of sorted survival times}
  \item{even.sort}{the event vector permuted using the permutation required to sort the survival times}
  \item{xval.surv}{n by n matrix of xvalidated survival curves, with rows indexed by time}
}
\references{

Aitkin, M. and Clayton, C. (1980), The Fitting of Exponential, Weibull and
Extreme Value Distributions to Complex Censored Survival Data using GLIM.
Applied Statistics, 29: 156-163.

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

Griffin, J. E., and Brown, P. J. (2005). Alternative prior distributions
for variable selection with very many more variables than observations.
Technical report available at
http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic/griffin/personal/

Cawley, G., and Talbot, N.,(2006)
Gene selection in cancer classification using sparse logistic
regression with Bayesian regularization
Bioinformatics 22,19,2 }

\author{Harri Kiiveri based on code by Glen Stone }
\seealso{HGcoxreg, HGsurv }
\examples{
set.seed(123)
# generate data
x <- matrix(rnorm(50 * 200), nrow = 50)
lp <- 3 + x[, 1]

# weibull shape parameter 1 - exponential distribution
a <- 1

# simulate survival times
time <- rweibull(50, a, exp(-lp / a))

# no censoring
event <- rep(1, 50)
# get cross validated survival curves
res<-xval.survival(x,time,event,fold=10,trace=TRUE)
# plot cross validated survival curves
matplot(res$times,res$xval.surv)

}
\keyword{ models }
\keyword{ survival }
\keyword{ cross validation}
