\name{plot.xval.HGpath}
\alias{plot.xval.HGpath}
\title{Plot HGpath cross validation error rates}
\description{
Plots the cross validation error over the calculated solution path
}
\usage{
\method{plot}{xval.HGpath}(obj, par = "delta", logpar = T,...)
}
\arguments{
  \item{obj}{ a xval.HGpath object }
  \item{par}{ one of "bbess" or "delta" defining which variable to plot on the axis }
  \item{logpar}{logical varible, TRUE means plot variable on log scale. }
  \item{\dots}{ any other parameters used by the method}
}
\details{
If one of length(pathk) or length(pathb) is one plot cross validated error rate and
standard deviation with "par" on x axis if length(pathk)=1. If both length(pathb) and length(pathk) are
greater than one produce a 2D contour plot
}
\value{
 Simply produces a plot
}
\author{Harri Kiiveri }
\seealso{HGpath, xval.HGpath, HGgaussian, HGmultc }
\examples{
#    res<-xval.HGpath(x,y,method=HGmultc)
#    plot(res)
#    plot(res,par="bbess",logpar=TRUE)
#    res$par.vals  # parameter combinations available for paths
}
\keyword{models }
\keyword{multivariate }
