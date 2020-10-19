AFFYMedianNorm <- function(X) {
	mds <- apply(X, 1, median)
	return(X / mds)
}

AFFYQuantileNorm <- function(X) {
	info <- 0

	temp <- .C("quantile_norm",
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		X = as.double(X),
		info = as.integer(info),
		PACKAGE = "RChip")

	if (info != 0)
		stop("Some problem in AFFYQuantileNorm")

	Xt <- matrix(temp$X, ncol = ncol(X), nrow = nrow(X))
	dimnames(Xt) <- dimnames(X)

	return(Xt)
}
