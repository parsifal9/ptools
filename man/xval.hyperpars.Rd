\name{xval.hyperpars}
\alias{xval.hyperpars}
\title{Cross validate the selection of hyperparameters}
\description{
  This function does a full cross validation of a model i.e it also cross validates
  the cross validated choice of the hyperparameters.
}
\usage{
xval.hyperpars(x, y, method = HGmultc, fold = NULL, foldbk = fold,
    trace = TRUE, etol=.05, bbess = 1e+07, kbess = seq(0, 1, 0.1),
    b0sc = 5, ...)
}
\arguments{
  \item{x}{ real matrix with n rows (samples) and p columns (variables or genes }
  \item{y}{ real n by 1 vector of response values appropriate to the method}
  \item{method}{ one of \code{\link{HGmultc}}, \code{\link{HGordcat}},
	\code{\link{HGgaussian}}, \code{\link{HGglm}}, \code{\link{HGcoxreg}},
	\code{\link{HGsurv}} }
  \item{fold}{ an integer specifying the number of folds to use in the cross validation.
	Default value NULL gives leave one out cross validation}
  \item{foldbk}{ an integer specifying the number of folds to use in the cross validation
                 for each combination of hyperparameters. Default value NULL gives
                 leave one out cross validation}
  \item{trace}{ a logical variable with value TRUE if information about cross validation
	progress to be output and value FALSE if no output required}
	\item{etol}{find models whose cross validation error is less than or equal to (1+etol) times the
              minimum cross validation error over all hyperparameter combinations. Must be >=0.
              See details below.}
  \item{bbess}{ vector of test values of scale parameter in sparsity prior,
	default 1e7.}
  \item{kbess}{ vector of test values of shape parameter in sparsitiy prior - 0 to 1,
	default 0 }
  \item{b0sc}{ vector of test values of length of initial beta vector - likely useful
	values in the range 1 - 25. Ignored in the methods \code{\link{HGmultc}},
	\code{\link{HGordcat}}}
  \item{\dots}{ any other arguments required by the chosen method }
}

\details{
   For each fold, the observations left in are used to find optimal values of the
   parameters bbess, kbess and b0sc. This is done by computing the foldbk cross validation error
   for each combination of bess, kbess and b0sc specified and chosing the "best" combination of parameters
   as follows. First identify all hyperparameter combinations with crossvalidated error rates within
   (1+etol) of the minimum. From this set find the subset of models with smallest value of kbess. Within this new set
   find the set of models with smallest value of b0sc (if present). Finally from this last set chose the model with
   largest value of bbess. Note that if etol=0 then the hyperparameters corresponding to the model with smallest
   cross validation error are chosen.  This process can be thought of as producing the "smallest" model with cross validation
   within (1+etol) of the minimum cross validation error.

   The folds used in  this process are chosen at random for each combination. The selected model is then
   used to predict the values of the observations held out. NOTE: this process can take a very long time for
   models with large numbers of variables as the number of model fits required is
   fold*foldbk*length(bbess)*length(kbess)*length(b0sc). Some time may be saved by not varying bbess.
}

\value{
   \item{xvfv }{a vector of fully cross validated fitted values}
  \item{grp }{a vector identifying which fold each observation was in}
  \item{details}{a matrix identifying the optimal choice of the parameters for each fold}
}
\references{ Ambroise, C. and McLachlan, G. J. (2002) Selection bias in gene extraction
on the basis of microarray gene-expression data. PNAS 2002, 99, 6562-6566 }

\author{Harri Kiiveri }
\seealso{ optim.hyperpar, xvalidate, DeleteRepeat}

\examples{
# EXAMPLE - linear regression
# generate  data
x <- matrix(rnorm(200 * 200), ncol = 200)
y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)

# do full cross validation of selection of bbess = c(1, 10),
# kbess = (0, 0.5), b0sc = c(1, 5)
tmp <- xval.hyperpars(x, y, method = HGgaussian, bbess = c(1, 10),
	kbess = c(0, 0.5), b0sc = c(1, 5), fold = 2)

# look at results
plot(tmp$xvfv, y)
print(tmp$details)
}

\keyword{ models }
