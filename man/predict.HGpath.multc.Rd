\name{predict.HGpath.multc}
\alias{predict.HGpath.multc}
\title{ Produce predicted values from HGpath object  }
\description{
Produce predicted (fitted) values for the entire solution path or rid from a HGpath object
}
\usage{
\method{predict}{HGpath.multc}(res, x)
}
\arguments{
  \item{res}{ a HGpath object }
  \item{x}{ the x matrix to use for prediction}
}
\details{
 Produces the predicted class for each model in the path or grid.
}
\value{
An n by length(pathb)*length(pathk) matrix of predicted values.
}
\author{Harri Kiiveri }
\seealso{ HGpath, HGgaussian, HGmultc  }
\examples{
# res<-HGpath(x,y,method=HGmultc,pathk=seq(1,0,-.01),pathb=c(.05,.05,1))
 # fv<-predict(res,x)
 # table(y,fv[,1])
}
\keyword{ models }
\keyword{ multivariate }
