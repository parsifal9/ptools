\name{fastRP}
\alias{fastRP}
\alias{rpartcallback}
\title{
Fast Recursive Partitioning and Regression Trees
}
\description{
Fit a \code{rpart} model
}
\usage{
fastRP(X, Y, weights = rep(1, nrow(X)), offset = NULL,
    RPmethod = "class", x = FALSE, y = TRUE, parms, control, cost,
    \dots)
}
\arguments{
\item{X}{ n by p matrix of expression values}
\item{Y}{ n by 1 vector of response values }
\item{weights}{ optional case weights }
\item{offset}{ask Glenn}
\item{RPmethod}{
one of \code{"anova"}, \code{"class"}
It is wisest to specifiy the method directly, especially as
more criteria are added to the function.

Alternatively, \code{method} can be a list of functions named
\code{init}, \code{split} and \code{eval}.
}
\item{x}{
keep a copy of the \code{x} matrix in the result.
}
\item{y}{
keep a copy of the dependent variable in the result. If missing and
\code{model} is supplied this defaults to \code{FALSE}.
}
\item{parms}{
optional parameters for the splitting function.
Anova splitting has no parameters.
For classification splitting, the list can contain any of:
the vector of prior probabilities (component \code{prior}), the loss matrix
(component \code{loss}) or the splitting index (component \code{split}).  The
priors must be positive and sum to 1.  The loss matrix must have zeros
on the diagnoal and positive off-diagonal elements.  The splitting
index can be \code{gini} or \code{information}.  The default priors are
proportional to the data counts, the losses default to 1,
and the split defaults to \code{gini}.
}
\item{control}{
options that control details of the \code{rpart} algorithm.
}
\item{cost}{
a vector of non-negative costs, one for each variable in the
model. Defaults to one for all variables.  These are scalings to be
applied when considering splits, so the improvement on splitting on a
variable is divided by its cost in deciding which split to choose.
}
\item{\dots}{
arguments to \code{rpart.control} may also be specified in the call to
\code{rpart}.  They are checked against the list of valid arguments. In particular
minsplit and minbucket are useful to set to prevent strange results
when sample sizes are small}
}

\value{
an object of class \code{rpart}, a superset of class \code{tree}.
}

\details{
This differs from the \code{tree} function mainly in its handling of surrogate
variables. In most details it follows Breiman et.\ al.\ quite closely. Glenn
has hacked this so that the model formulae and associated overhead have been removed.
It should work fine with large numbers of variables.
}

\references{
Breiman, Friedman, Olshen, and Stone. (1984)
\emph{Classification and Regression Trees.}
Wadsworth.
}

\seealso{
\code{\link{rpart.control}}, \code{\link{rpart.object}},
\code{\link{summary.rpart}}, \code{\link{print.rpart}}
}

\examples{
# classification example
w <- abs(sin(1:200))
x <- matrix(rnorm(200 * 200), nrow = 200)
y <- (x[, 1] > 0) + 1

fit <- fastRP(x, y, RPmethod = "class")
fv <- predict(fit, x, type = "class")
table(y, fv)

# regression exxample
x <- matrix(rnorm(200 * 200), ncol = 200)
y <- 2 * x[, 1] + 1 + 0.1 * rnorm(200)

fit1 <- fastRP(x, y, RPmethod = "anova")
fv1 <- predict(fit1, x, type = "anova")

plot(y, fv1)

# otherwise on some devices the text is clipped
par(mfrow = c(1, 1), xpd = NA)
plot(fit)
text(fit, use.n = TRUE)
plot(fit1)
text(fit1, use.n = TRUE)
}
\keyword{tree}

