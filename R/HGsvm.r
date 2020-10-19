HGsvm <- function (x, f, weights = rep(1, nrow(x)), sparsity.prior = "NG",
	bbess = 1e7, kbess = 0,  scale = 10, initb = NULL, no.prior = 1)
{
   scale<--abs(scale) # dont update scale - used for hinge approx
   # force response to be -1 and 1
   a<-sort(unique(f))
   if(length(a)!=2 | ((a[1]!=1)&(a[2]!=2))){
   cat("error: response must have values 1 or 2","\n")
   return()
   }
   f[f==a[1]]<--1
   f[f==a[2]]<-1
	# add column of ones and get initial value
	x <- cbind(rep(1, nrow(x)), x)
	if (is.null(initb[1])) {
		bhat <- HGinit.g(x, f, eps = 0, lb = 1, weights, b0sc)
		nl <- max(abs(bhat))
		bhat <- bhat/nl
	}
	else {
		bhat <- initb
  	}

	#
	# call engine
	#
	res <- HGengine(x, f , event = NULL, weights = weights,
		bbess = bbess, kbess = kbess, scale = scale, bhat = bhat,
		Le = Le.svm, DLde = DLde.svm, D2Lde = D2Lde.svm, HGsc = HGsc.svm,
		sparsity.prior = sparsity.prior,
		no.prior = no.prior, control = HGcontrol(tolc = 1e-3))

	# post process output                          +
   #sigma2 <- res$scale/ysc^2
	# if (scale < 0)
	#	sigma2 <- abs(scale)
	#
	#compute alpha
	tmp<-lm(res$beta*bbess~t(x)-1)
  alpha<-tmp$coefficients*f
	zz <- list(beta = res$beta, S = res$S, fv = res$eta,
		varids = which(res$S)-1, sigma2 = scale, alpha=alpha)
	class(zz) <- c("HGsvm", "HG")
	zz
}

#-functions required by HGengine for svm 2 class model

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

Le.svm <- function(eta,y,event,weights,scale)
{
 lp<-(1-y*eta)*scale
 lp<-pmax(pmin(lp,200),-200)
 L<- -sum(log(1+exp(lp))*weights)/scale
 L
}

DLde.svm <- function(eta,y,event,weights,scale)
{
 lp<-(1-y*eta)*scale
 lp<-pmax(pmin(lp,200),-200)
 p<-exp(lp)/(1+exp(lp))
 y*p*weights
}

D2Lde.svm <- function(eta,y,event,weights,scale)
{
 lp<-(1-y*eta)*scale
 lp<-pmax(pmin(lp,200),-200)
 p<-exp(lp)/(1+exp(lp))
 -scale*p*(1-p)*weights
}

HGsc.svm<- function(eta,y,event,weights,scale)
{
  scale
}
