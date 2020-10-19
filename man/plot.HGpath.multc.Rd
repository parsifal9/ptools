\name{plot.HGpath.multc}
\alias{plot.HGpath.multc}
\title{ HGmultc solution path plot}
\description{
Plots the solution path for a HGpath object with method HGmultc
}
\usage{
\method{plot}{HGpath.multc}(obj, class = 1, cols = NULL, par = "delta", logpar = T,...)
}

\arguments{
  \item{obj}{ a HGpath object }
  \item{class}{ which class (1,...,G) parameters to plot path for }
  \item{cols}{a vector of rows of obj$par.vals to use to define the path. If NULL
  plot default paths as defined in the details below.   }
  \item{par}{one of "delta" or "bbess" defining x axis variable  }
  \item{logpar}{logical variable, if TRUE plot log(parameter) on x axis }
  \item{\dots}{ any other parameters used by the method}
}
\details{
If cols=NULL and either one of length(bpath) or length(kpath) is one, plots
 the solution path as a function of the bpath or kpath variable which is not of length one.
 If length(pathb) and length(kpath) are both greater than one and cols=NULL plot the diagonal path
 indexed by (bbess[1:q],kbess[1:q]) where q is min(length(pathb),length(pathk)). Otherwise use
 the path sequence specified by cols
}
\value{
  Produces a plot of the solution path
}
\author{Harri Kiiveri }
\seealso{HGpath, xval.HGpath and HGmultc }
\examples{
#    res<-HGpath(x,y,method=HGmultc)
#    plot(res)
#    plot(res,class=2,par="bbess")
#    res$par.vals  # parameter combinations available for paths
}
\keyword{ models }
\keyword{ multivariate }
