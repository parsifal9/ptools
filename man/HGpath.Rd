\name{HGpath}
\alias{HGpath}
\title{ GeneRave solution paths}
\description{
Computes solution paths for HGgaussian and HGmultc over a grid of values of bbess and kbess
}
\usage{
HGpath(x, y, method, weights=rep(1,nrow(x)), pathk = seq(1, 0, -0.01), pathb = c(0.1, 1e-04, 50), trace = TRUE, control = HGcontrol(), ...)
}
\arguments{
  \item{x}{ an n by p data matrix }
  \item{y}{ a n by 1 response vector }
  \item{method}{one of HGgaussian or HGmultc }
  \item{weights}{an n by 1 vector of observation weights}
  \item{pathk}{ a decreasing sequence of kbess values starting from one and going to zero }
  \item{pathb}{ a vector of three components, the first the largest values of bbess on a path,
  the second the smallest value of bbess, and the third the total number of bbess values which will be
   equally spaced on a log scale  }
  \item{trace}{ if TRUE report progress over path values}
  \item{control}{ control parameters for HGgaussian and HGmultc}
  \item{\dots}{ any other relevant arguments for HGmultc or HGgaussian }
}
\details{
  Uses the fast algorithm of Friedman, Hastie and Tibshirani(2008) to compute the initial
  solution path for an L1 fit. Warm starts are then used by HGgaussian and HGmultc to efficiently
  compute all the parameter estimates over the grid of bbess and kbess values.
}
\value{
A list with components
  \item{bpath }{for method=HGgaussian a p+1 by length(pathb)*length(pathk) sparse matrix of parameter estimates (including the constant)
  for the grid or path. For method=HGmultc a p+1 by length(pathb)*length(pathk)*G sparse matrix of parameter
  estimates for the grid or path, where G is the number of classes. The column indexing is such that in order of
  increasing speed of variation we have bbess, kbess and class }
  \item{par.vals }{a (length(pathb)*length(pathk)) by 3 matrix of parameter values defining the path or grid.
  The first column gives the bbess values, the second column the kbess values and the third column the
  delta=(2/bbess)^0.5 values. For L1 regression delta is n*lambda. }
}
\references{ Friedman,J., Hastie,T., and Tibshirani, R.(2008) Regularization paths for generalized linear
models via coordinate descent. Manuscript submitted for publication. }
\author{Harri Kiiveri }
\note{ There are predict and plot methods for the object produced by this function.
}
\seealso{ HGgaussian and HGmultc }
\examples{
# Compute L1 constrained regression solution path (kbess fixed at 1)
  x<-matrix(rnorm(5000),nrow=50,ncol=100)
  lp<-x[,1]+2*x[,3]+6*x[,4]
  y<-lp+rnorm(50)*.1
  res<-HGpath(x,y,trace=FALSE,method=HGgaussian,pathk=1)
  plot(res)
  res<-HGpath(x,y,trace=FALSE,method=HGgaussian,pathk=1, pathb=c(.1,1e-6,50))
  plot(res,par="bbess",logpar=TRUE)
# get numbers of variables in each model
  nv<-colSums(res$bpath!=0)
# Compute solution path for classification problem  (bbess fixed at a value)
  phat<-exp(lp)/(1+exp(lp))
  yind<-rbinom(50,1,phat)+1
  res<-HGpath(x,yind,trace=FALSE,method=HGmultc,pathk=seq(1,0,-.01), pathb=c(1e-1,1e-1,1))
  plot(res)
}
\keyword{models }
\keyword{multivariate }
