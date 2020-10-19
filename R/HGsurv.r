HGsurv <- function (x, time, event, weights = rep(1, nrow(x)),
	sparsity.prior = "NG", bbess = 1e7, kbess = 0, b0sc = 1, scale = -1,
	initb = "FALSE", loginthaz = NULL, aupdate = NULL, a = NULL, no.prior = 1)
{
# new version of code
# loginthaz is function for computing log integrated hazard function
# aupdate - function to update estmates of parameters in log integrated
# hazard function
# a -intial values for log integrated hazard function parameter(s)
# these must all be specified to run this function

	if (is.null(loginthaz) || is.null(aupdate) || is.null(a))
		stop("Error: you must specify all three of loginthaz, aupdate and a.")

	lb <- 1

	# add column of ones and get initial value
	x <- cbind(rep(1, nrow(x)), x)

	if (initb[1] == FALSE) {
		#logt <- log(time)
		bhat <- HGinit.s(x, time, event, lb, weights, b0sc, a, ofn = loginthaz)
		nl <- sqrt(sum(bhat * bhat))

		if (b0sc < 0)
			b0sc <- nl

		ysc <- b0sc / nl
		bhat <- bhat * ysc
	}
	else {
		bhat <- initb
		bmax <- max(abs(x %*% bhat))

		# some protection against inappropriate intial values
		if (bmax > 20)
			bhat <- bhat * 20 / bmax

		#nl <- sqrt(sum(bhat * bhat))
		#ysc <- b0sc / nl
		#bhat <- bhat * ysc
		#scale <- scale * ysc^2
	}
 	#
	# call engine - note that event is the response in the poisson model!!!
	#
	res <- HGengine(x, event, event = time, weights = weights, bbess = bbess,
		kbess = kbess, scale = scale, bhat = bhat,
		Le = Le.p, DLde = DLde.p, D2Lde = D2Lde.p, HGsc = HGsc.p,
		sparsity.prior = sparsity.prior, no.prior = no.prior,
		Offset = loginthaz, Offset.par.update = aupdate, a = a)

	# post process output
	fhat <- exp(res$eta) > 0.5
	lp <- x %*% res$beta
	lfv <- loginthaz(event, time, res$a)
	surv0 <- exp(-exp(lfv))
	survall <- outer(exp(lfv), as.vector(exp(lp)), FUN = "*")
	survall <- exp(-survall)
	surv0 <- cbind(time, surv0, event)
	colnames(surv0) <- c("event time", "baseline survival fn", "censor vector")
	survall <- cbind(time, survall)

	cnames <- paste("subject", 1:nrow(x))
	colnames(survall) <- c("event time", cnames)

	k <- order(time)

	zz <- list(beta = res$beta, S = res$S, lp = res$eta, surv0 = surv0[k,], survall = survall[k,],
		varids = which(res$S) - 1, haz.par = res$a, fv = exp(res$eta))
	class(zz) <- c("HGsurv", "HG")
	zz
}

# additional functions required by HGsurv

HGinit.s <- function (x, time, event, lb, weights, b0sc, a, ofn)
{
	y <- as.vector((event * 100 - (1 - event) * 100) - 100 - ofn(event, time, a))
	x <- sweep(x, 1, sqrt(weights), FUN = "*")
	y <- y * sqrt(weights)

	if (nrow(x) < ncol(x)) {
		xtx <- x %*% t(x)
		yy <- xtx %*% y
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		z <- backsolve(qr.R(temp), rhs, nrow(xtx))
		bss <- (t(x) %*% (y - z)) / lb
	}
	else {
		xtx <- t(x) %*% x
		yy <- t(x) %*% y
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		bss <- backsolve(qr.R(temp), rhs, nrow(xtx))
	}
	bss
}
HGinit.s.arch <- function (x, time, event, lb, weights, b0sc, a, ofn)
{
# old version for historical reasons no longer used
	y <- as.vector(log(event * 0.99 + (1 - event) * 0.01) - ofn(event, time, a))
	x <- sweep(x, 1, sqrt(weights), FUN = "*")
	y <- y * sqrt(weights)

	if (nrow(x) < ncol(x)) {
		xtx <- x %*% t(x)
		yy <- xtx %*% y
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		z <- backsolve(qr.R(temp), rhs, nrow(xtx))
		bss <- (t(x) %*% (y - z)) / lb
	}
	else {
		xtx <- t(x) %*% x
		yy <- t(x) %*% y
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		bss <- backsolve(qr.R(temp), rhs, nrow(xtx))
	}

	bss
}

Le.p <- function(eta, y, event, weights, scale)
{
	sum(weights * (y * eta - exp(eta)))
}

DLde.p <- function(eta, y, event, weights, scale)
{
	weights * (y - exp(eta))
	# grad(Le.p,eta,y=y,event=event,weights=weights,scale=scale)-test code
}

D2Lde.p <- function(eta, y, event, weights, scale)
{
	-weights * exp(eta)
	# diag(hessian(Le.p,eta,y=y,event=event,weights=weights,scale=scale)) - test code
}

HGsc.p <- function(eta, y, event, weights, scale)
{
	scale <- 1
	scale
}

# weibull functions - for testing
wloginthaz <- function (event, time, a)
{
	a * log(time)
}

waupdate <- function (event, time, a0, eta, weights)
{
	p0 <- exp(eta)
	a1 <- sum(log(time) * (p0 - event) * weights)
	a1 <- sum(weights * event) / a1
	a0 + 0.1 * (a1 - a0)
}


