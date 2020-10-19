\name{HGsvm}
\alias{HGsvm}
\title{HG support vector machine classification}
\description{
Fits a support vector machine using a differentiable approximation to the hinge loss formulation
and a sparsity prior
}
\usage{
HGsvm(x, f, weights = rep(1, nrow(x)), sparsity.prior = "NG", bbess = 1e+07, kbess = 0, scale = 10, initb = NULL, no.prior = 1)
}
\arguments{
  \item{x}{ real n by p data matrix of predictors}
  \item{f}{ real n by 1 vector of classes labelled 1 or 2 }
  \item{weights}{ real n by 1 vector of observation weights }
  \item{sparsity.prior}{one of "NG", "NEG", "DEG", "SCAD", or "R" }
  \item{bbess}{ real number in the range [0,1e7]. When kbess equals one for the NG prior
  the usual penalty parameter in L1 regression is (2/bbess)^0.5. In this case values
  of bbess in the interval (0,1] will often be useful }
  \item{kbess}{real number in the range [0,1]. The value for the SCAD prior can
  be bigger than one, see for example the help for HGmultc.}
  \item{scale}{ parameter determining the accuracy of the  approximation to the hinge loss criterion}
  \item{initb}{if NULL generate inital value automatically, if set equal to a p+1 by 1
		vector use this vector to initialise iterations. }
  \item{no.prior}{Default is 1, don't penalise the intercept. }
}
\details{
Fits a support vector machine using the approximation |1-x|+ =((scale)^-1)*log(1+exp(scale*(1-x))).
Larger values of scale make the approximation closer. The default should be adequate for most cases.
}
\value{
An object of class HGsvm which is a list with the following components:
\item{beta}{ p+1 by 1 vector of parameter estimates   }
\item{S}{
	p+1 by 1 logical vector with value true if variable selected i.e
	the corresponding beta is nonzero
}
\item{fv}{n by 1 vector of fitted class values  }
\item{varids}{identifiers of variables selected, 0 denotes intercept }
\item{sigma2}{value of parameter in hingle loss approximation  }
\item{alpha}{usual alpha parameters in SVM, see reference}
}
\references{Zhang,T., and Frank J. Oles, F. J.(2001) Text categorization based on
regularized linear classification methods. Information Retrieval,4, 1, 5-31.

Hastie,T., Rosset,S., Tibshirani,R., Zhu,J. (2004) The entire
regularization path for the support vector machine.
Journal of Machine Learning Research, 5, 1391-1415.}
\author{ Harri Kiiveri }
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
#generate data
x<-matrix(rnorm(5000),nrow=50,ncol=100)
  lp<-1+x[,1]+2*x[,3]+6*x[,4]
   y<-lp+rnorm(50)*.1
  phat<-exp(lp)/(1+exp(lp))
  yind<-rbinom(50,1,phat)+1
# fit model
  res<-HGsvm(x,yind,kbess=0.5,bbess=0.1)
  fv<-predict(res,x)
  table(yind,fv)
# do tenfold cross validation - not run
#  tmp<-xvalidate(x,yind,method=HGsvm,kbess=0.5,bbess=0.1,fold=10,trace=T)
# table(yind,tmp)
}
\keyword{models}
\keyword{multivariate}
