\name{predict.HGpath.gaussian}
\alias{predict.HGpath.gaussian}
\title{ Produce predicted values from HGpath object}
\description{
Produce predicted (fitted) values for the entire solution path or rid from a HGpath object
}
\usage{
\method{predict}{HGpath.gaussian}(res, x)
}
\arguments{
  \item{res}{ a HGpath object }
  \item{x}{ the x matrix to use for prediction }
}
\details{
Produces the linear predictor for each model on the path or grid.
}
\value{
An n by length(pathb)*length(pathk) matrix of predicted values.
}
\author{Harri Kiiveri }
\seealso{ HGpath, HGgaussian, HGmultc }
\examples{
 # res<-HGpath(x,y,method=HGgaussian,pathk=seq(1,0,-.01),pathb=c(.05,.05,1))
 # fv<-predict(res,x)
 # plot(y,fv[,1])
}
\keyword{ models }
\keyword{ multivariate }
