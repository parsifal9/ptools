dlda <- function(X, y, priors) {
	n <- nrow(X)
	p <- ncol(X)
	G <- max(y)

	pp <- .C("dldawrap",
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(X),
		as.integer(y),
		means = matrix(0.0, G, p),
		vars = double(p),
		counts = integer(G),
		DUP = FALSE,
		PACKAGE = "RChip"
	)

	if (missing(priors)) {
		priors <- rep(1.0,G)
		priors[which.max(pp$counts)] <- 1.0001
		priors <- priors / sum(priors)
	}

	pp <- c(pp[c("means", "vars", "counts")], list(priors = priors))

	class(pp) <- "dlda"
	return(pp)
}

dqda <- function(X, y, priors) {
	n <- nrow(X)
	p <- ncol(X)
	G <- max(y)

	pp <- .C("dqdawrap",
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(X),
		as.integer(y),
		means = matrix(0.0, G, p),
		vars = matrix(0.0, G, p),
		counts = integer(G),
		DUP = FALSE,
		PACKAGE = "RChip")

	if (missing(priors))   {
		priors <- rep(1.0, G)
		priors[which.max(pp$counts)] <- 1.0001
		priors <- priors / sum(priors)
	}

	pp <- c(pp[c("means", "vars", "counts")], list(priors = priors))

	class(pp) <- "dqda"
	return(pp)
}

predict.dlda <- function(object, newdata, type = c("class", "prob"), ...)
{
	type <- match.arg(type)
	G <- length(object$counts)
	n <- nrow(newdata)
	p <- ncol(object$means)

	pp <- .C("preddldawrap",
		as.integer(p),
		as.integer(G),
		as.double(object$means),
		as.double(object$vars),
		as.integer(object$counts),
		as.double(object$priors),
		as.integer(n),
		as.double(newdata),
		R = matrix(0.0, n, G),
		DUP = FALSE,
		PACKAGE = "RChip"
	)

	if (type == "prob")
		return(pp$R)
	else
		return(apply(pp$R, 1, which.max))
}

predict.dqda <- function(object, newdata, type = c("class", "prob"), ...)
{
	type <- match.arg(type)
	G <- length(object$counts)
	n <- nrow(newdata)
	p <- ncol(object$means)

	pp <- .C("preddldawrap",
		as.integer(p),
		as.integer(G),
		as.double(object$means),
		as.double(object$vars),
		as.integer(object$counts),
		as.double(object$priors),
		as.integer(n),
		as.double(newdata),
		R = matrix(0.0, n, G),
		DUP = FALSE,
		PACKAGE = "RChip"
	)

	if (type == "prob")
		return(pp$R)
	else
		return(apply(pp$R, 1, which.max))
}

sdda <- function(obj, ...) {
	if (is.matrix(obj))
		return(sdda.matrix(obj, ...))
	else
		return(sdda.dda(obj, ...))
}

sdda.matrix <- function(obj, y, priors, ...) {
	if (missing(priors))
		tt <- dlda(obj, y)
	else
		tt <- dlda(obj, y, priors)

	return(sdda.dda(tt, obj, y, ...))
}

sdda.dda <- function(obj, X, y, start = rep(FALSE, p), never = rep(FALSE, p),
	useprob = TRUE, usecache = TRUE, usexval = TRUE, jmax = -20, ...)
{
	p <- ncol(X)
	n <- nrow(X)
	G <- length(obj$counts)

	if (inherits(obj,"dlda"))
		lda <- 1
	else
		lda <- 0

	pflg <- ifelse(useprob,1,0)
	cflg <- ifelse(usecache,1,0)
	xflg <- ifelse(usexval,1,0)

	if (any(start))
		never[start] <- FALSE

	jj <- max(1, abs(jmax))

	pp <- .C("sddawrap",
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(X),
		as.integer(y),
		as.double(obj$means),
		as.double(obj$vars),
		as.integer(obj$counts),
		as.double(obj$priors),
		as.integer(start),
		as.integer(never),
		S = logical(p),
		as.integer(lda),
		as.integer(pflg),
		as.integer(cflg),
		as.integer(xflg),
		as.integer(jmax),
		ecrit = double(jj),
		pcrit = double(jj),
		DUP = FALSE,
		PACKAGE = "RChip"
	)

	obj$S <- pp$S
	names(obj$S) <- colnames(X)
	obj$ecrit <- pp$ecrit
	obj$pcrit <- pp$pcrit

	class(obj) <- c("sdda", class(obj))
	return(obj)
}

predict.sddanull <- function(object, newdata, type = c("class", "prob"), ...)
{
	type <- match.arg(type)
	G <- length(object$counts)
	priors <- object$priors

	n <- nrow(newdata)
	R <- matrix(priors, nrow = n, ncol = G, byrow = TRUE)

	if (type == "prob")
		return(R)
	else
		return(apply(R, 1, which.max))
}

predict.sdda <- function(object, newdata, ...) {
	S <- object$S
	if (all(!S)) return(predict.sddanull(object, newdata, ...))

	objT <- list(
		means = object$means[ , S, drop = FALSE],
		vars =
			if (inherits(object,"dlda")) object$vars[S]
			else object$vars[ , S, drop = FALSE],
		counts = object$counts,
		priors = object$priors
	)

	class(objT) <- class(object)[-1]
	predict(objT, newdata[ , S, drop = FALSE], ...)
}

plotdiag <- function(obj, ...) UseMethod("plotdiag")

plotdiag.sdda <- function(obj, ...) {
	k <- min(sum(obj$S) + 1, length(obj$ecrit))
	op <- par(mfrow = c(2, 1))
	plot(1:k, obj$pcrit[1:k])
	plot(1:k, obj$ecrit[1:k])
	par(op)
	invisible()
}

plot.sdda <- function(x, X, y, ...) {
	S <- x$S
	j <- sum(S)

	if (j == 0) {
		warning("No genes selected - cannot plot")
		return(invisible())
	}

	if (j == 1) {
		plot(as.factor(y), X[, S])
	} else {
		require(MASS)
		ll <- lda(X[, S], y)
		xx <- X[,S] %*% ll$scaling
		PairsPlot(xx, lab = y, col = y + 1)
	}

	invisible()
}

which.genes.sdda <- function(object, ...) return(which(object$S))
