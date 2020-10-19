HGengine <- function(X, y, event, weights = rep(1, nrow(X)), bbess = 1e+07, kbess = 0,
	scale, bhat, Le, DLde, D2Lde, HGsc, sparsity.prior = "NG", no.prior = NULL,
	Offset = NULL, Offset.par.update = NULL, a = NULL, control = HGcontrol())
{
	# X - n by p data matrix
	# y - n by 1 response vector
	# event - working vector - required for survival analysis (set to NULL usually)
	# weights - n by 1 vector of positive weights
	# bbess, kbess- hyperparameters in prior
	# scale - scale parameter - negative if fixed
	# bhat -inital estimate of beta
	# Le(eta,y,event,weights,scale) - function to compute likelihood as a function of eta=X'beta
	# DLde - function to compute first derivative of likelihood as function of eta
	# D2Lde - function to compute second derivative of likelihood as fn of eta
	# or expected second derivative or diagonal of derivative or some negative (semi) definite matrix
	# can be a n by 1 vector (diags) or an n by n neg def matrix
	# HGcexp - function to compute conditional expectation
	# no.prior - column ids of variables which are not constrained by the prior
	# HGsc - function for updating scale parameter - fn of (eta,y,event,cwts,scale)
	# Offset - function to compute vector to add to linear predictor
	# Offset.par.update - function to update offset parameters during outer iterations
	# a - initial values of parameters a in offset function
	# control - a list of control parameters, see HGcontrol
	scale.fix <- FALSE
	if (scale < 0) {
		scale.fix <- TRUE
		scale <- abs(scale)
	}

	tolc <- control$tolc

	HGcexp <- match.prior(sparsity.prior)

	n <- nrow(X)
	p <- ncol(X)

	weights <- abs(weights)
	nn <- sum(weights > 0)
	weights <- (weights/sum(weights)) * nn

	# compute intial offset if present
	offset <- NULL
	if (!is.null(Offset))
		offset <- Offset( y, event, a)

	S <- bhat != 0
	S[no.prior] <- TRUE
	R <- HGeta(X, bhat, S, offset, control$eta.lim)
	Q1 <- HGpost(y, event, bhat, bhat, S, R, weights, bbess, kbess, Le, scale,
		HGcexp = HGcexp, no.prior)
	beta <- bhat
	j <- 1

	repeat {
		Q0 <- Q1
		repeat {
			delta <- HGdirection(X, y, event, beta, bhat,
				S, R, weights, bbess, kbess, DLde, D2Lde, scale,
				HGcexp = HGcexp, no.prior)
			alpha <- 1
			nc <- 0
			repeat {
				btmp <- beta
				tt <- S
				btmp[tt] <- btmp[tt] + alpha * delta

				R <- HGeta(X, btmp, S, offset, control$eta.lim)
				Q2 <- HGpost(y, event, btmp, bhat, S, R, weights, bbess,
					kbess, Le, scale, HGcexp = HGcexp, no.prior)
				#cat(alpha,Q2-Q1,sum(S),(Q2>Q1-tolc),"\n")
				if (Q2 > Q1 - tolc)
					break
				nc <- nc + 1
				if (nc > 100)
					stop("Error in line search. No upward step could be found")
				alpha <- alpha/2
			}
			#cat(nc,alpha,Q2,Q1,"\n")
			beta <- btmp

			if (abs(Q2 - Q1) < tolc)
				break

			Q1 <- Q2

			if (alpha * max(abs(delta)) < control$tolerance) {
				break
			}
		}
		j <- j + 1
		if (exists("Qlog"))
			Qlog <<- c(Qlog, Q1)
		S <- (abs(beta) > control$epsilon)
		S[no.prior] <- TRUE
		bhat <- beta

		if (abs(Q1 - Q0) < control$tolerance)
			break

		#R <- HGeta(X, bhat, S, offset)

		if (!scale.fix)
			scale <- HGsc(R, y, event, weights, scale)   # moved from after line below

			# update offset if present
		if (!is.null(Offset)) {
			a <- Offset.par.update(y, event, a, R, weights)
			offset <- Offset(y, event, a = a)
		}
			R <- HGeta(X, bhat, S, offset, control$eta.lim)   #added

		Q1 <- HGpost(y, event, bhat, bhat, S, R, weights, bbess, kbess,
			Le, scale, HGcexp = HGcexp, no.prior)



#		cat("outer", j, sum(S), "\n") #- test

		if (j > control$maxit)
			break
	}

	zz <- list(beta = bhat, S = S, eta = R, scale = scale, a = a)
	class(zz) <- c("HGengine", "HG")
	zz
}

HGdirection <- function (X, y, event, beta, bhat, S, R, weights, bbess, kbess, DLde, D2Lde, scale,
	HGcexp, no.prior)
{
	tt <- S
	bssi <- HGcexp(bbess, kbess, bhat[tt])
	bssi <- bssi^-0.5
	pg <- sum(tt)

	if (pg == 0)
		return(NULL)

	n <- nrow(X)

	if (!is.null(no.prior)) {
		kappa <- 10
		k <- which(tt)
		ll <- match(no.prior, k)
		bssi[ll] <- kappa
		beta[tt][ll] <- 0
	}

	B <- X[, tt, drop = FALSE] * matrix(bssi, nrow = n, ncol = pg, byrow = TRUE)
	am <- D2Lde(R, y, event, weights, scale)
	#if(is.matrix(am)) A <- chol(-am)%*%B

	if (is.square.matrix(am)) {
		#cat("here I am","\n") - test
		svdam <- try(svd(-am), silent = TRUE)

		if (class(svdam) == "try-error") {
			svdam <- eigen(-am)
			svdam$d <- svdam[[1]]
			svdam$u <- svdam[[2]]
		}

		D <- svdam$d
		ind <- (D > 1e-08)
		U <- svdam$u[, ind]
		D <- sqrt(D[ind])
		A <- crossprod(U, B)
		A <- sweep(A, 1, D, FUN = "*")
	}
	else {
		A <- as.vector(sqrt(-D2Lde(R, y, event, weights, scale))) * B
	}

	B <- crossprod(B, DLde(R, y, event, weights, scale)) - beta[tt] / bssi

	if (pg < n) {
		A <- crossprod(A)
		diag(A) <- diag(A) + 1
		Z <- solve(A, B)
		delta <- bssi * Z
	}
	else {
		M <- crossprod(t(A))
		diag(M) <- diag(M) + 1
		Z <- solve(M, A %*% B)
		delta <- bssi * (B - crossprod(A, Z))
	}
	#cat(delta[1:5],"\n")
	delta
}

HGcexpold.ng <- function (b, k, beta)
{
	gi <- beta^2/(2 * b)
	y <- 2 * sqrt(gi)
	e <- beta * beta
	e <- 1/e
	e <- e * y * besselK(y, 1.5 - k)/besselK(y, 0.5 - k)
	e
}

HGcexp.ng<-function (b, k, beta)
{
    delta <- sqrt(2/b)
    if(k==0){
    ab<-abs(beta)
    e<-1/ab^2 +delta/ab}
    else if (k == 1) e <- delta/abs(beta)
    else {
        y <- delta * abs(beta)
        e <- beta * beta
        e <- 1/e
        z<- besselK(y, 1.5 - k)/besselK(y, 0.5 - k)
        z[is.na(z)]<-1  # fix for underflow
        e <- e * y *z
    }
    e
}

HGpost <- function (y, event, beta, bhat, S, R, weights, bbess, kbess, Le, scale,
	HGcexp, no.prior)
{
	Qf <- Le(R, y, event, weights, scale)
	S[no.prior] <- FALSE
	bssi <- HGcexp(bbess, kbess, bhat[S])
	bssi <- bssi^-0.5
	Qf <- Qf - 0.5 * sum((beta[S]/bssi)^2, na.rm = TRUE)
	if (!is.numeric(Qf))
		Qf <- -1e32
	Qf
}

HGeta <- function (X, beta, S, offset, eta.lim)
{
	R <- X[ , S, drop = FALSE] %*% beta[S]

	if (!is.null(offset))
		R <- R + offset

	cupcap(R, eta.lim)
}

cupcap <- function (x, lim)
{
	pmin(pmax(x, -lim), lim)
}

is.square.matrix <- function(a)
{
	ind <- FALSE
	a <- as.matrix(a)
	if (nrow(a) == ncol(a))
		ind <- TRUE
	ind
}
