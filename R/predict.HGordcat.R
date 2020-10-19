"predict.HGordcat" <- function(object, newdata, type = "class", ...)
{
	#browser()
	G <- ncol(object$P)
	y <- rep(1, nrow(newdata))

	res <- orderedcat(newdata, y, G, weights = y)
	beta <- c(object$theta, object$beta)

	lp <- res$xs %*% beta
	lp <- pmax(-100, lp)
	lp <- pmin(100, lp)
	lp <- as.matrix(lp)

	p <- exp(lp)
	p <- p / (1 + p)
	res <- oddstop(p, G)

	if (type == "class") {
		f <- res$lab
	}
	else {
		f <- res$p
	}

	f
}
