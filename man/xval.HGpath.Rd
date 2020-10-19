\name{xval.HGpath}
\alias{xval.HGpath}
\title{ Solution path cross validation}
\description{
Cross validates the solution paths produced by HGpath
}
\usage{
xval.HGpath(x, y, method = HGgaussian, pathk = seq(1, 0, -0.01), pathb = c(0.1, 1e-04, 50), fold = 10, trace = T, 
control = HGcontrol(tolerance = 1e-04, tolc = 0.01),weights=rep(1,nrow(x)), ...)
}
\arguments{
  \item{x}{an n by p data matrix   }
  \item{y}{a n by 1 response vector  }
  \item{method}{one of HGgaussian or HGmultc  }
  \item{pathk}{a decreasing sequence of kbess values starting from one and going to zero  }
  \item{pathb}{a vector of three components, the first the largest values of bbess on a path,
  the second the smallest value of bbess, and the third the total number of bbess values which will be
   equally spaced on a log scale   }
  \item{fold}{number of folds in the cross validation  }
  \item{trace}{ if TRUE report progress over path values }
  \item{control}{control parameters for HGgaussian and HGmultc  }
  \item{weights}{a vector of observation weights}
  \item{\dots}{ any other relevant arguments for HGmultc or HGgaussian}
}
\details{
 Computes cross validated fitted values and cross validated error rates and associated standard deviations
 for the grid or path specified by pathb and pathk
}
\value{
A list with components
  \item{xvfv }{a n by (length(pathb)*length(pathk)) matrix of cross validated fitted values}
  \item{par.vals }{a (length(pathb)*length(pathk)) by 3 matrix of parameter values defining the path or grid.
  The first column gives the bbess values, the second column the kbess values and the third column the
  delta=(2/bbess)^0.5 values. For L1 regression delta is n*lambda.}
  \item{grp}{a n by 1 vector identifying the folds used in the cross validation}
  \item{xve}{a length(pathb)*length(pathk) vector of cross validated error rates}
  \item{xvsd}{a length(pathb)*length(pathk) vector of cross validated standard deviations}
  \item{jmin}{a vector identifying the rows of par.vals above which have cross validated errorr rates less
  than or equal to the minimum cross validated error rate plus one standard deviation}
}
\author{Harri Kiiveri }
\note{ There is a plot method for the object produced by this function
}
\seealso{ HGpath, HGgaussian and HGmultc }
\examples{
 # Cross validate L1 regression solution path (kbess fixed at 1)
  x<-matrix(rnorm(5000),nrow=50,ncol=100)
  lp<-x[,1]+2*x[,3]+6*x[,4]
  y<-lp+rnorm(50)*.1
  res<-xval.HGpath(x,y,trace=FALSE,method=HGgaussian,pathk=1, pathb=c(1,1e-6,25))
  plot(res)
  # plot cross validated fitted values versus observed for one model
  k<-res$jmin[1]
  plot(res[[1]][,k],y)
}
\keyword{models }
\keyword{multivariate}
