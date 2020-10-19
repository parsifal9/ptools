HGordcat <- function (x, y, weights = rep(1, nrow(x)), sparsity.prior = "NG",
	bbess = 1e7, kbess = 0)
{
	no.prior <- 1:(max(y) - 1)
	lb <- 1
	scale <- 1
	# set maximum iterations
	p <- ncol(x)

	G <- max(y)
	if (G < 1) {
		stop("Must be more than 2 classes.")
	}

	weights <- abs(weights)
	nn <- sum(weights > 0)
	weights <- (weights/sum(weights)) * nn

	# preprocess data
	res <- HGmake.data.oc(x, y, G, weights)

	# compute intial value
	e <- matrix(nrow = length(res$ys), ncol = 2)
	for (i in 1:2) {
		e[, i] <- as.numeric(res$ys + 1 == i)
	}

	bhat <- HGinit.oc(res$xs, res$y, lb, weights)

	# call engine
	res1 <- HGengine(res$xs, res$ys, event = NULL, weights = res$weights,
		bbess = bbess, kbess = kbess, scale = scale, bhat = bhat,
		Le = Le.l, DLde = DLde.l, D2Lde = D2Lde.l, HGsc = HGsc.l,
		sparsity.prior = sparsity.prior, no.prior = no.prior)

	# post process data
	expeta <- exp(res1$eta)
	pd <- 1+expeta
	p <- cbind(expeta/pd,1/pd)
	out <- HGoddstop(p, G, col = 2)

	i <- G:1
	q <- apply(out$p[, i], 1, cumsum)
	q <- out$p/t(q[i, ])

	# unpack results
	theta <- res1$beta[1:(G - 1)]
	K <- G:length(res1$beta)
	coeffs <- res1$beta[K]
	S <- res1$S[K]

	zz <- list(beta = -coeffs, S = S, P = out$p, varids = which(S),
		theta = -theta, class = out$lab)
	class(zz) <- c("HGordcat", "HG")
	zz
}

HGinit.oc <- function (x, e, lb, weights)
{
	y <- log(e * 18 + (1 - e)/18 )
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
	nl <- sum(bss^2)
	bss <- bss * 15/nl
}

Le.l <- function(eta, y, event, weights, scale)
{
	sum(weights*(y*eta-log(1+exp(eta))))
}

DLde.l <- function(eta, y, event, weights, scale)
{
	#a <- grad(Le.l,eta,y=y,event=event,weights=weights,scale=scale)
	expeta <- exp(eta)
	weights * (y - expeta / (1 + expeta))
}


D2Lde.l <- function(eta, y, event, weights, scale)
{
	#b <- hessian(Le.l,eta,y=y,event=event,weights=weights,scale=scale)
	expeta <- exp(eta)
	p0 <- 1 / (1 + expeta)
	p <- expeta * p0
	-weights * p0 * p
}

HGsc.l <- function(eta, y, event, weights, scale)
{
	scale <- 1
	scale
}

HGmake.data.oc <- function (x, y, G, weights)
{
	n <- nrow(x)
	G1 <- G - 1
	tmp <- matrix(rep(weights, G1), nrow = n)
	weights <- as.vector(t(tmp))
	xs <- matrix(0, nrow = n * G1, ncol = G1 + ncol(x))
	I <- diag(G1)
	i <- rep(1:n, rep(G1, n))
	xs[, G:ncol(xs)] <- x[i, ]
	for (j in 1:G1) {
		xs[, j] <- rep(I[, j], n)
	}
	I <- matrix(0, nrow = n, ncol = G)
	for (j in 1:n) {
		I[j, y[j]] <- 1
	}
	I <- apply(I, 1, cumsum)
	ys <- as.vector(I[1:G1, ])
	list(xs = xs, ys = ys, weights = weights)
}

HGoddstop <- function (p, G, col = 1)
{
	G <- G - 1
	n <- nrow(p)
	m <- n/G
	pc <- matrix(0, nrow = m, ncol = (G + 1))
	a <- matrix(0, nrow = m, ncol = G + 1)

	for (i in 1:G) {
		a[, i + 1] <- p[seq(i, n, G), col]
	}

	pc[, G + 1] <- a[, G + 1]

	if (G >= 2) {
		for (k in (G + 1):2) {
			pc[, k - 1] <- (a[, k - 1]/a[, k]) * pc[, k] * (1 -a[, k])
		}
	}

	pc[, 1] <- 1 - pc %*% rep(1, (G + 1))
	lab <- apply(pc, 1, which.max)

	list(p = pc, lab = lab)
}
