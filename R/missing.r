completeGene <- function(Xmat) {
	tt <- apply(Xmat, 2, function(x) any(is.na(x)))
	return(Xmat[, !tt])
}

completeArray <- function(Xmat) {
	tt <- apply(Xmat, 1, function(x) any(is.na(x)))
	return(Xmat[!tt, ])
}

impute <- function(Xmat) {
	rmeans <- rowMeans(Xmat, na.rm = TRUE)
	cmeans <- colMeans(Xmat, na.rm = TRUE)
	gmean <- mean(Xmat, na.rm = TRUE)

	isNA <- is.na(Xmat)

	n <- nrow(Xmat)
	p <- ncol(Xmat)
	impval <- matrix(rmeans, ncol = p, nrow = n, byrow = FALSE) +
		matrix(cmeans, ncol = p, nrow = n, byrow = TRUE) - gmean

	Xmat[isNA] <- impval[isNA]
	Xmat
}
