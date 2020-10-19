HGglm <- function (x, f, event = NULL, weights = rep(1, nrow(x)),
	sparsity.prior = "NG", bbess = 1e7, kbess = 0, b0sc = 15, scale = -1,
	initb = "FALSE", model = "N",
	fvalfn = NULL, ifvalfn = NULL, varfn = NULL, drfn = NULL, devfn = NULL,
	scale.updatefn = NULL, no.prior = 1)
{
	# model = N - normal disribution
	#     P - poisson distribution
	#     B - binomial - logistic regression
	#     G - gamma
	#     IG - inverse gaussian distribution
	#     Own - define own model (you must specify your own intial value in initb)
	# fvalfn - a function of the linear predictor eta which returns mu - the fitted value
	# ifvalfn - function of the fitted value which returns the linear predictor (inverse of fvalfn)
	# varfn - a function of mu which returns the value of the variance function
	# drfn - function of mu which returns dmu/deta
	# devfn - a function of mu, y, event, scale, weights which returns the value of the deviance
	# scale.updatefn - a function of mu, y, event, scale, weights which returns the value of the scale parameter

	# add column of ones and get initial value
	x <- cbind(rep(1, nrow(x)), x)

	#
	# fit normal distribution --------------------------------------------------
	if (model == "N") {
	#
		if (initb[1] == FALSE) {
			bhat <- HGinit.g(x, f, eps = 0, lb = 1, weights, b0sc)
			nl <- sqrt(sum(bhat * bhat))

			if (b0sc < 0)
				b0sc <- nl

			ysc <- b0sc / nl
			bhat <- bhat * ysc
		}
		else {
			bhat <- initb
			#nl <- sqrt(sum(bhat * bhat))
			#ysc <- b0sc / nl
			#bhat <- bhat * ysc
			#scale <- scale * ysc^2
			ysc <- 1
		}

		#
		# call engine
		#
		res <- HGengine(x, f * ysc, event = NULL, weights = weights,
			bbess = bbess, kbess = kbess, scale = scale, bhat = bhat,
			Le = Le.g, DLde = DLde.g, D2Lde = D2Lde.g, HGsc = HGsc.g,
			sparsity.prior = sparsity.prior, no.prior = no.prior)

		# post process output
		sigma2 <- res$scale / ysc^2

		if (scale < 0)
			sigma2 <- abs(scale)

		zz <- list(beta = res$beta / ysc, S = res$S, fv = res$eta / ysc,
			varids = which(res$S) - 1, scale = sigma2, model = model)
	}

	#
	# fit binomial distribution ------------------------------------------------
	#
	if (model == "B") {
		if (is.null(event))
			event <- rep(1, length(f))

		if (any(f > event)) {
			stop("Error: one or more y values exceeds the number of binomial trials.")
		}

		if (initb[1] == FALSE) {
			bhat <- HGinit.b(x, f, event, lb = 1, weights)
		}
		else {
			bhat <- initb
		}

		#
		# call engine
		#
		# binomial sample sizes in event, default to 1

		res <- HGengine(x, f, event, weights = weights, bbess = bbess,
			kbess = kbess, scale = scale, bhat = bhat,
			Le = Le.b, DLde = DLde.b, D2Lde = D2Lde.b, HGsc = HGsc.b,
			sparsity.prior = sparsity.prior, no.prior = no.prior)

		expeta <- exp(res$eta)
#		p0 <- event * expeta / (1 + expeta)
		p0 <- expeta / (1 + expeta)

		# post process output
		zz <- list(beta = res$beta, S = res$S, fv = p0,
			varids = which(res$S) - 1, scale = res$scale, model = model,
			bd = event)
	}

	#
	# fit poisson distribution ------------------------------------------------
	#
	if (model == "P") {
		if (initb[1] == FALSE) {
			bhat <- HGinit.p(x, log(f + .01), lb = 1, weights)
			nl <- sqrt(sum(bhat * bhat))

			if (b0sc < 0)
				b0sc <- nl

			ysc <- b0sc / nl
			bhat <- bhat * ysc
		}
		else {
			bhat <- initb
		}

		#
		# call engine
		#
		res <- HGengine(x, f, event, weights = weights, bbess = bbess,
			kbess = kbess, scale = scale, bhat = bhat,
			Le = Le.p, DLde = DLde.p, D2Lde = D2Lde.p, HGsc = HGsc.p,
			sparsity.prior = sparsity.prior, no.prior = no.prior)

		# post process output
		zz <- list(beta = res$beta, S = res$S, fv = exp(res$eta),
			varids = which(res$S) - 1, scale = scale, model = model)
	}

	#
	# fit gamma distribution
	#
	if (model == "G") {
		if (initb[1] == FALSE) {
			bhat <- HGinit.p(x, 1 / (f + .01), lb = 1, weights)
			nl <- sqrt(sum(bhat * bhat))

			if (b0sc < 0)
				b0sc <- nl

			ysc <- b0sc / nl
			bhat <- bhat * ysc
		}
		else {
			bhat <- initb
		}

		if (any(x %*% bhat <= 0)) {
			stop("Error in initial value: linear predictor has negative values.")
		}

		#
		# call engine
		#
		res <- HGengine(x, f, event, weights = weights, bbess = bbess,
			kbess = kbess, scale = scale, bhat = bhat,
			Le = Le.ga, DLde = DLde.ga, D2Lde = D2Lde.ga, HGsc = HGsc.ga,
			sparsity.prior = sparsity.prior, no.prior = no.prior)

		# post process output
		zz <- list(beta = res$beta, S = res$S, fv = 1 / res$eta,
			varids = which(res$S) - 1, scale = scale, model = model)
	}

	#
	# fit inverse gaussian distribution
	#
	if (model == "IG") {
		if (initb[1] == FALSE) {
			bhat <- HGinit.p(x, 1 / (f + .01)^2, lb = 1, weights)
			nl <- sqrt(sum(bhat * bhat))

			if (b0sc < 0)
				b0sc <- nl

			ysc <- b0sc / nl
			bhat <- bhat * ysc
		}
		else {
			bhat <- initb
		}

		if (any(x %*% bhat <= 0)) {
			stop("Error in initial value: linear predictor has negative values")
		}

		#
		# call engine
		#
		res <- HGengine(x, f, event, weights = weights, bbess = bbess,
			kbess = kbess, scale = scale, bhat = bhat,
			Le = Le.ig, DLde = DLde.ig, D2Lde = D2Lde.ig, HGsc = HGsc.ig,
			sparsity.prior = sparsity.prior, no.prior = no.prior)

		# post process output
		zz <- list(beta = res$beta, S = res$S, fv = 1 / (res$eta)^2,
			varids = which(res$S)- 1, scale = scale, model = model)
	}

	if (model == "Own") {
		#
		# fit own model-------------------------------------------------------------
		#

		# Set up environments so they can access the relevant user-defined
		# functions
		environment(Le.o)   <- environment()
		environment(DLde.o)  <- environment()
		environment(D2Lde.o) <- environment()
		environment(HGsc.o)  <- environment()
		#environment(fvalfn)<-environment()  # add other own functions here ?

		if (initb[1] == FALSE) {
		bhat <- HGinit.g(x, ifvalfn(f), eps = 0, lb = 1, weights, b0sc)
			#stop("You must supply an initialisation vector when using your own model.")
			nl <- sqrt(sum(bhat * bhat))

			if (b0sc < 0)
				b0sc <- nl

			ysc <- b0sc / nl
			bhat <- bhat * ysc
		}
		else {
			bhat <- initb
		}

		#
		# call engine
		#
		res <- HGengine(x, f, event = NULL, weights = weights, bbess = bbess,
			kbess = kbess, scale = scale, bhat = bhat,
			Le = Le.o, DLde = DLde.o, D2Lde = D2Lde.o, HGsc = HGsc.o,
			sparsity.prior = sparsity.prior, no.prior = no.prior)

		# post process output
		zz <- list(beta = res$beta, S = res$S, fv = fvalfn(res$eta),
			varids = which(res$S) - 1, scale = res$scale, model = model)
	}

	class(zz) <- c("HGglm", "HG")
	return(zz)
}

# other functions required for HGglm

Le.o <- function(eta, y, event, weights, scale)
{
	#browser()
	mu <- fvalfn(eta)
	devfn(mu, y, event, scale, weights)
}

DLde.o <- function(eta, y, event, weights, scale)
{
	#browser()
	#a <- grad(Le.o,eta,y=y,event=event,weights=weights,scale=scale)
	mu <- fvalfn(eta)
	(weights / scale) * (drfn(mu) / varfn(mu)) * (y - mu)
}

D2Lde.o <- function(eta, y, event, weights, scale)
{
	#browser()
	mu <- fvalfn(eta)
	-(weights / scale) * (drfn(mu)^2 / varfn(mu))
}

HGsc.o <- function(eta, y, event, weights, scale)
{
	mu <- fvalfn(eta)
	scale.updatefn(mu, y, event, weights, scale)
}

# binomial functions
Le.b <- function(eta, y, event, weights, scale)
{
	sum(weights * (y * eta - event * log(1 + exp(eta))))
}

DLde.b <- function(eta, y, event, weights, scale)
{
	#a <- grad(Le.l,eta,y=y,event=event,weights=weights,scale=scale)
	expeta <- exp(eta)
	weights * (y - event * expeta / (1 + expeta))
}

D2Lde.b <- function(eta, y, event, weights, scale)
{
	#b <- hessian(Le.l,eta,y=y,event=event,weights=weights,scale=scale)
	expeta <- exp(eta)
	p0 <- 1 / (1 + expeta)
	p <- expeta * p0
	-weights * event * p0 * p
}

HGsc.b <- function(eta, y, event, weights, scale)
{
	scale <- 1
	scale
}

HGinit.b <- function (x, e, event, lb, weights)
{
#	if (is.null(event))
#		event <- 1

	y <- e / event
	y <- pmax(y, .001)
	y <- pmin(y, .999)
	y <- log(y / (1 - y))
	ws <- sqrt(weights)
	x <- sweep(x, 1, ws, FUN = "*")
	y <- y * ws
	if (nrow(x) < ncol(x)) {
		xtx <- x %*% t(x)
		yy <- xtx %*% y
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		z <- backsolve(qr.R(temp), rhs, nrow(xtx))
		bss <- (t(x) %*% (y - z))/lb
	}
	else {
		xtx <- t(x) %*% x
		yy <- t(x) %*% y
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		bss <- backsolve(qr.R(temp), rhs, nrow(xtx))
	}
	#bss <- sweep(bss, 2, apply(abs(bss), 2, max), "/")
	#nl <- sum(bss^2)
	#bss <- bss * 15 / nl
	bss
}

# poisson functions

HGinit.p <- function (x, f, lb, weights)
{
	# eps not required at
	x <- sweep(x, 1, sqrt(weights), FUN = "*")
	y <- f * sqrt(weights)

	if (nrow(x) < ncol(x)) {
		xtx <- crossprod(t(x))
		yy <- crossprod(xtx, y)
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		z <- backsolve(qr.R(temp), rhs, nrow(xtx))
		bss <- (t(x) %*% (y - z)) / lb
	}
	else {
		xtx <- crossprod(x)
		yy <- crossprod(x, y)
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
#	grad(Le.p,eta,y=y,event=event,weights=weights,scale=scale)-test code
}

D2Lde.p <- function(eta, y, event, weights, scale)
{
	-weights * exp(eta)
#	diag(hessian(Le.p,eta,y=y,event=event,weights=weights,scale=scale)) - test code
}

HGsc.p <- function(eta, y, event, weights, scale)
{
	scale <- 1
	scale
}

# functions for gamma distribution
Le.ga <- function(eta, y, event, weights, scale)
{
	if (any(eta <= 0))
		L <- -1e32
	else
		L <- scale * sum(weights * (log(eta) - eta * y))

	L
}

DLde.ga <- function(eta, y, event, weights, scale)
{
	eta <- pmax(1e-10, eta)
	scale * weights * (1 / eta - y)
#	grad(Le.p,eta,y=y,event=event,weights=weights,scale=scale)-test code
}

D2Lde.ga <- function(eta, y, event, weights, scale)
{
	eta <- pmax(1e-10, eta)
	-scale * weights / eta^2
#	diag(hessian(Le.p,eta,y=y,event=event,weights=weights,scale=scale)) - test code
}

HGsc.ga <- function(eta, y, event, weights, scale)
{
	scale <- sum((y * eta - 1)^2) / length(y)
	scale <- 1 / scale
}

# functions for inverse gaussian distribution
Le.ig <- function(eta, y, event, weights, scale)
{
	if (any(eta <= 0))
		L <- -1e32
	else
		L <- sum(weights * y * (sqrt(eta) - 1 / y))

	-L
}

DLde.ig <- function(eta, y, event, weights, scale)
{
	one <- rep(1, length(y))
	eta <- pmax(1e-10, eta)
	-(weights * y * (one - 1 / (y * sqrt(eta))))
#	grad(Le.p,eta,y=y,event=event,weights=weights,scale=scale)-test code
}

D2Lde.ig <- function(eta, y, event, weights, scale)
{
	eta <- pmax(1e-10, eta)
	-weights / (2 * eta^1.5)
#	diag(hessian(Le.p,eta,y=y,event=event,weights=weights,scale=scale)) - test code
}

HGsc.ig <- function(eta, y, event, weights, sc)
{
	fv <- 1 / sqrt(eta)
	linv <- sum(wts * ((y - fv) / (y * fv))^2 * y) / sum(weights)
	#cat(linv,"\n")
	linv <- max(linv, .01)
	1 / linv
	# sc<-sc+0.5*(sc1-sc)
}

# functions to test "own" models ---------------
# gaussian fns
#fvalfn <- function(eta) {eta}

#varfn <- function(mu) {1}

#drfn <- function(mu) {1}

#devfn <- function(mu, f, event, scale, weights)
#{
#	# f is y variable
#	# p0 is the linear predictor at point b0
#	#
#	L <- -0.5 * sum(weights) * log(scale) -
#		0.5 * sum(((weights * (f - mu)^2) / scale))
#	L
#}

#scupdate <- function(mu, f, event, scale, weights)
#{
#	scale <- sum(weights * (f - mu)^2) / length(f)
#	scale <- max(scale, 1e-2)
#	scale
#}

# logistic fns

#fvalfn <- function (eta)
#{
#	eta <- cupcap(eta, 100)
#	exp(eta) / (1 + exp(eta))
#}

#devfn <- function (mu, y, event, scale, weights)
#{
#	#browser()
#	mu <- pmin(1 - 1e-10, mu)
#	mu <- pmax(1e-10, mu)
#	sum(weights * y * log(mu) + weights * (1 - y) * log(1 - mu))
#}

#drfn <- function (mu)
#{
#	mu <- pmin(1 - 1e-10, mu)
#	mu <- pmax(1e-10, mu)
#	mu * (1 - mu)
#}

#varfn <- function (mu)
#{
#	mu <- pmin(1 - 1e-10, mu)
#	mu <- pmax(1e-10, mu)
#	mu * (1 - mu)
#}

#scupdate <- function (mu, y, event, scale, weights)
#{
#	scale
#}
