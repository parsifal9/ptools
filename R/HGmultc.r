HGmultc<-function (X, y, weights = rep(1, nrow(X)), initb=NULL, sparsity.prior = "NG",
    lambda = 1, bbess = 1e+07, kbess = 0, adj.prior = TRUE, control = HGcontrol())
{
    if (max(y) == 1) {
        cat("Error: class labels must be 1,2,...,G", "\n")
        return()
    }
    tol0 <- control$tolerance
    tol1 <- tol0
    eps0 <- control$epsilon
    eps1 <- 1e-4 # changed from eps0
    tolc <- control$tolc
    HGIcexpb <- match.prior(sparsity.prior)
    n <- nrow(X)
    p <- ncol(X)
    X <- cbind(rep(1, n), X)
    p <- p + 1
    G <- max(y)
    G1 <- G
    ff <- 1
    if (adj.prior)
        ff <- (G - 1)/G
    weights <- abs(weights)
    nn <- sum(weights > 0)
    weights <- (weights/sum(weights)) * nn
    if(is.null(initb))bhat <- HGIinitialisew(X, y, eps1, lambda, weights)
    else bhat<-initb
    S <- bhat != 0
    S[1,]<-TRUE # fix so constant always present even if not in b0
    R <- HGIprob(X, bhat, S)
    Q1 <- HGIQfunwb(y, bhat, bhat, S, R, weights, bbess, kbess,
        ff, HGIcexpb)
    beta <- bhat
    j <- 1
    repeat {
        Q0 <- Q1
        repeat {
            delta <- vector("list", G1)
            for (g in 1:G1) {
                delta[[g]] <- HGIdirectionwb(X, y, beta, bhat,
                  S, R, g, weights, bbess, kbess, ff, HGIcexpb)
            }
            alpha <- 1
            nc <- 0
            repeat {
                btmp <- beta
                for (g in 1:G1) {
                  tt <- S[, g]
                  btmp[tt, g] <- btmp[tt, g] + alpha * delta[[g]]
                }
                R <- HGIprob(X, btmp, S)
                Q2 <- HGIQfunwb(y, btmp, bhat, S, R, weights,
                  bbess, kbess, ff, HGIcexpb)
                if (Q2 > Q1 - tolc)
                  break
                nc <- nc + 1
                #cat(Q1,Q2,nc,tolc,"\n")
                if (nc > 100)
                  stop("Error in line search\n No upward step could be found")
                alpha <- alpha/2
            }
            beta <- btmp
            if (abs(Q2 - Q1) < tolc)
                break
            Q1 <- Q2
            if (alpha * max(abs(unlist(delta))) < tol1)
                break
        }
        j <- j + 1
        if (exists("Qlog"))
            Qlog <<- c(Qlog, Q1)
        S <- (abs(beta) > eps0)
        S[1, ] <- TRUE
        bhat <- beta
        if (abs(Q1 - Q0) < tol0)
            break
        R <- HGIprob(X, bhat, S)
        Q1 <- HGIQfunwb(y, bhat, bhat, S, R, weights, bbess,
            kbess, ff, HGIcexpb)
        if (j > control$maxit) {
            cat("Iteration Limit exceeded.\n")
            break()
        }
    }
    #class <- apply(R, 1, which.max)
    class<-max.col(R)
    zz <- list(beta = bhat, S = S, P = R, class = class)
    class(zz) <- c("HGmultc", "HG")
    zz
}

HGIQfunwb.old <- function (y, beta, bhat, S, R, weights, bbess, kbess, ff,
	HGIcexpb)
{
	n <- nrow(R)
	G <- max(y)
	#Rt <- 1 - apply(R, 1, sum)
	Rt<-1-rowSums(R)
	Qf <- sum(weights * log(pmax(1e-100, R[cbind(1:n, y)])))
	S[1, ] <- FALSE

	for (g in 1:ncol(S)) {
		bssi <- HGIcexpb(bbess, kbess, bhat[S[, g], g])
		bssi <- bssi^-0.5
		Qf <- Qf - ff * 0.5 * sum((beta[S[, g], g] / bssi)^2, na.rm = TRUE)#added ff
	}

	Qf
}

HGIQfunwb <- function (y, beta, bhat, S, R, weights, bbess, kbess, ff,
	HGIcexpb)
{
	n <- nrow(R)
	G <- max(y)
	#Rt <- 1 - apply(R, 1, sum)
	Rt<-1-rowSums(R)
	Qf <- sum(weights * log(pmax(1e-100, R[cbind(1:n, y)])))
	S[1, ] <- FALSE
  #
  bhat<-as.vector(bhat*S)  # need S here for equivalence
  beta<-as.vector(beta)
  ind<-(bhat!=0)
  bhat<-bhat[ind]
  bssi<-HGIcexpb(bbess,kbess,bhat)
  bssi<-bssi^-0.5
  Qf<-Qf-ff*0.5*sum((beta[ind]/bssi)^2,na.rm=T)
	Qf
}

HGIdirectionwb <- function (X, y, beta, bhat, S, R, g, weights, bbess, kbess,
	ff, HGIcexpb)
{
	tt <- S[, g]

	bssi <- HGIcexpb(bbess, kbess, bhat[tt, g])
	bssi <- bssi^-0.5

	pg <- sum(tt)
	if (pg == 0)
		return(NULL)

	n <- nrow(X)
	kappa <- 10
	bssi[1] <- kappa
	beta[tt, g][1] <- 0
	B <- X[, tt, drop = FALSE] * matrix(bssi, nrow = n,
		ncol = pg, byrow = TRUE) / ff^0.5 #added ff^0.5 here
	A <- sqrt(weights * R[, g] * (1 - R[, g])) * B
	B <- crossprod(B, weights * ((y == g) - R[, g])) - ff * beta[tt, g] / bssi #added ff

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

	delta
}

HGIinitialisew <- function (X, y, eps, lambda, weights)
{
	n <- nrow(X)
	p <- ncol(X)
	G <- max(y)
	R <- matrix(log(eps), nrow = n, ncol = G)
	R[cbind(1:n, y)[y <= G, ]] <- -log(eps)

	R <- sweep(R, 1, weights, FUN = "*")
	X <- sweep(X, 1, sqrt(weights), FUN = "*")

	if (p <= n) {
		A <- crossprod(X)
		diag(A) <- diag(A) + lambda
		B <- crossprod(X, R)
		beta <- solve(A, B)
	}
	else {
		A <- crossprod(t(X))
		B <- A %*% R
		diag(A) <- diag(A) + lambda
		beta <- crossprod(X, R - solve(A, B))/lambda
	}

	beta <- t(scale(t(beta), center = TRUE, scale = FALSE))
	beta <- sweep(beta, 2, apply(abs(beta), 2, max), "/")
	25 * beta
}

#HGIcexpb <- function(b, k, beta)
#{
#b scale parameter
# k shape parameter
# beta conditioning values for E{v-2 | beta}
# returns conditional expected value for expanded hg model
#
#	gi <- beta^2 / (2 * b)
#	y <- 2 * sqrt(gi)
#	e <- beta * beta
#	e <- 1 / e
#	e <- e * y * besselK(y, 1.5 - k) / besselK(y, 0.5 - k)
#	e
#}

HGIprob.old <- function (X, beta, S)
{
	R <- matrix(0, nrow = nrow(X), ncol = ncol(beta))

	for (g in 1:ncol(beta)) {
		tt <- S[, g]
		R[, g] <- exp(cupcap(X[, tt, drop = FALSE] %*%
			beta[tt, g, drop = FALSE], 100))
	}

	RG <- 1 / (apply(R, 1, sum))
	R <- sweep(R, 1, RG, "*")
	R
}

HGIprob <- function (X, beta, S)
{
#browser()
	R <- matrix(0, nrow = nrow(X), ncol = ncol(beta))

  R<-exp(cupcap(X%*%(beta*S),100))
  #R<-exp(cupcap(tcrossprod(X,t(beta*S)),100))

	#RG <- 1 / (apply(R, 1, sum))
	RG<-1/rowSums(R)
	#R <- sweep(R, 1, RG, "*")
	R<-R*RG
	R
}
