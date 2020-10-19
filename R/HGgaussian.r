HGgaussian<-function (x, f, weights = rep(1, nrow(x)), sparsity.prior = "NG",
    bbess = 1e+07, kbess = 0, b0sc = 5, scale = -1, initb = NULL,
    no.prior = 1, control = HGcontrol(eta.lim = 1e30))
{
    lambda<-1
    if(nrow(x)>ncol(x)){
    lambda<-1e-4
    b0sc<- -1    # approximate regression estimates in this case
    }
    x <- cbind(rep(1, nrow(x)), x)
    if (is.null(initb[1])) {
        bhat <- HGinit.g(x, f, eps = 0, lb = lambda, weights, b0sc)
        nl<-max(abs(bhat[-1]))
        if (b0sc < 0)
            b0sc <- nl
        ysc <- b0sc/nl
        bhat <- bhat * ysc
    }
    else {
        bhat <- initb
        ysc <- 1
    }
    res <- HGengine(x, f * ysc, event = NULL, weights = weights,
        bbess = bbess, kbess = kbess, scale = scale, bhat = bhat,
        Le = Le.g, DLde = DLde.g, D2Lde = D2Lde.g, HGsc = HGsc.g,
        sparsity.prior = sparsity.prior, no.prior = no.prior,
        control = control)
    sigma2 <- res$scale/ysc^2
    if (scale < 0)
        sigma2 <- abs(scale)
    zz <- list(beta = res$beta/ysc, S = res$S, fv = res$eta/ysc,
        varids = which(res$S) - 1, sigma2 = sigma2)
    class(zz) <- c("HGgaussian", "HG")
    zz
}

#-functions required by HGengine for regression model

HGinit.g <- function (x, f, eps, lb, weights, b0sc)
{
#  eps not required at all ?? b0sc
	x <- sweep(x, 1, sqrt(weights), FUN = "*")
	y <- f * sqrt(weights)
	if (nrow(x) < ncol(x)) {
		xtx <- crossprod(t(x))
		yy <- crossprod(xtx, y)
		xtx <- xtx + diag(nrow(xtx)) * lb
		temp <- qr(xtx)
		rhs <- qr.qty(temp, yy)
		z <- backsolve(qr.R(temp), rhs, nrow(xtx))
		bss <- (t(x) %*% (y - z))/lb
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

Le.g <- function(eta,y,event,weights,scale)
{
	n <- sum(weights)
	L <- -n*log(scale)/2 -0.5*sum(weights*(y-eta)^2)/scale
	L
}

DLde.g <- function(eta,y,event,weights,scale)
{
	weights*(y-eta)/scale
}

D2Lde.g <- function(eta,y,event,weights,scale)
{
	-weights/scale
}

HGsc.g <- function(eta,y,event,weights,scale)
{
  #s0<-scale
	scale <- sum(weights*(y-eta)^2)/sum(weights)
	scale <- max(scale,1e-4)
	#s0+0.1*(scale-s0)
}
