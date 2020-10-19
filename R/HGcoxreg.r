HGcoxreg <- function (x, time, event, weights = rep(1, nrow(x)),
	sparsity.prior = "NG", bbess = 1e7, kbess = 0, b0sc = 1, scale = -1,
	initb = FALSE, no.prior = NULL)
{
	# new version
	lb <- 1

	# preprocess data
	n <- nrow(x)

	# break ties in survival time
	k <- order(time, decreasing = FALSE)
	d <- abs(diff(time[k]))
	sds <- sqrt(min(d[d > 0]))
	ts <- time[k]

	if (any(d == 0)) {
		i <- (2:n)[d == 0]
		ts[i] <- ts[i] + rnorm(length(i), 0, sds / 3)
		time[k] <- ts
	}

	# reorder data
	k <- order(time, decreasing = FALSE)
	x <- x[k,]
	time <- time[k]
	event <- event[k]

	if (initb[1] == FALSE) {
		logt <- log(time)
		bhat <- HGinit.cr(x, event, logt, lb, weights, b0sc)
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
		#ysc <- b0sc/nl
		#bhat <- bhat * ysc
		#scale <- scale * ysc^2
	}
  	# refine intial value using ridge likelihood -omit
		refine <- FALSE
		if (refine) {
			res <- HGengine(x, time, event, weights = weights,
				bbess = 1, kbess = kbess,
				scale = scale, bhat = bhat,
				Le = Le.cr, DLde = DLde.cr, D2Lde = D2Lde.cr, HGsc = HGsc.cr,
				sparsity.prior = "R", no.prior = no.prior,
				control = HGcontrol(maxit = 2))

			bhat <- res$beta
		}
	#
	# call engine
	#
	res <- HGengine(x, time, event, weights = weights, bbess = bbess,
		kbess = kbess, scale = scale, bhat = bhat,
		Le = Le.cr, DLde = DLde.cr, D2Lde = D2Lde.cr, HGsc = HGsc.cr,
		sparsity.prior = sparsity.prior, no.prior = no.prior)

	# post process output
	res1 <- HGsurv.fn(event, res$eta)
	fv0 <- cbind(time, res1[[1]], event)
	colnames(fv0) <- c("event time", "baseline survival fn", "censor vector")

	fvall <- cbind(time, t(res1[[2]]))
	cnames <- paste("subject ", 1:n, sep = "")
	colnames(fvall) <- c("event time", cnames)

	# undo order in output
	lp <- res$eta
	lp[k] <- res$eta
	stmp <- fvall
	stmp[, c(1, k + 1)] <- fvall

	z <- list(beta = res$beta, S = res$S, lp = lp, surv0 = fv0, survall = stmp,
		varids = which(res$S))
	class(z) <- c("HGcoxreg", "HG")
	z
}

#-functions required by HGengine for cox regression model
HGinit.cr <- function (x, f, logt, lb, weights, b0sc)
{
	y <- as.vector((f * 100- (1 - f) * 100) - logt - 100)
	x <- sweep(x, 1, sqrt(weights), FUN = "*")
	y <- y * sqrt(weights)

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
	bss
}

HGinit.cr.archive <- function (x, f, logt, lb, weights, b0sc)
{
# old version kept here for historical reasons
	y <- as.vector(log(f * 0.999 + (1 - f) * 0.001) - logt)
	x <- sweep(x, 1, sqrt(weights), FUN = "*")
	y <- y * sqrt(weights)

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

	bss
}

Le.cr <- function(eta, y, event, weights, scale)
{
	n <- length(eta)
	Sj <- cumsum(exp(eta)[n:1])[n:1]
	L <- sum(event*(eta-log(Sj)))
	L
}

#
DLde.cr <- function(eta, y, event, weights, scale)
{
	#browser()
	n <- length(eta)
	expeta <- exp(eta)
	Sj <- cumsum(exp(eta)[n:1])[n:1]
	Rj <- (event/Sj)
	dL <- event-cumsum(Rj)*expeta
	# grad(Le.cr,eta,y=y,event=event,weights=weights,scale=scale)-test code
	dL
}

#
D2Lde.cr <- function(eta, y, event, weights, scale)
{
	#browser()
	n <- length(eta)
	W <- matrix(1, nrow = n, ncol = n)
	W <- upper.tri(W, diag = TRUE)
	expeta <- exp(eta)
	Sj <- cumsum(exp(eta)[n:1])[n:1]
	W <- sweep(W, 1, Sj, FUN = "/")
	W <- sweep(W, 2, expeta, FUN = "*")
	# calculate diagonals to add
	w <- crossprod(W, event)
	# calculate first component of second derivative
	W <- sweep(W, 1, event, FUN = "*")
	D2L <- crossprod(W)
	diag(D2L) <- diag(D2L) - w
	D2L
	# diag(hessian(Le.cr,eta,y=y,event=event,weights=weights,scale=scale)) - test code
}

#
HGsc.cr <- function(eta, y, event, weights, scale)
{
	scale <- 1
	scale
}

HGsurv.fn <- function (event, eta)
{
	n <- length(eta)
	expeta <- exp(eta)
	Sj <- cumsum(exp(eta)[n:1])[n:1]
	EPSI <- (as.vector(c(rep(1, n)) - (matrix(expeta,ncol = 1) / (Sj))))^exp(-eta)
	hatBASELINE <- EPSI
	hatBASELINE[1] <- EPSI[1]^event[1]

	for (i in 1:(n - 1)) {
		hatBASELINE[i + 1] <- hatBASELINE[i] * EPSI[i + 1]^event[i + 1]
	}

	MM <- matrix(hatBASELINE, nrow = 1)
	MM <- t(matrix(rep(MM, length(eta)), nrow = length(eta)))
	indivsurv <- apply(MM, 2, "^", expeta)
	list(hatBASELINE = hatBASELINE, indivsurv = indivsurv)
}


