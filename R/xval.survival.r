xval.survival <- function (X, y, event, method = HGcoxreg, fold = NULL,
	trace = FALSE, weights = rep(1, nrow(X)), ...)
{
#	if (missing(breaks) || is.null(breaks)) {
#		stop("argument 'breaks' must be specified")
#	}

# break time ties outside this function first
	i <- order(y)
	d <- abs(diff(y[i]))
	if (any(d == 0)) { stop("Error: break ties in survival times first") }

	if (is.null(colnames(X)))
		colnames(X) <- 1:ncol(X)
	n <- nrow(X)
	C <- rep(0, n)

	xsurv <- matrix(NA, nrow = n, ncol = n)

	if (is.null(fold)) {
		grp <- 1:n
	}
	else {
		grp <- rep(1:fold, n/fold + 1)[1:n]
		grp[order(as.numeric(y) + runif(n))] <- grp
	}

	for (i in sort(unique(grp))) {
		if (trace)
			cat("X-validate group", i, "of", max(grp), "\n")

		tt <- (grp == i)

		if (is.null(event))
			tmp <- method(X[!tt, ], y[!tt], weights = weights[!tt], ...)
		else
			tmp <- method(X[!tt, ], y[!tt], event[!tt], weights = weights[!tt], ...)

		lp <- predict(tmp, X[tt, , drop = FALSE], type = type)
		C[tt] <- lp

# get cross validated survival curves
		nn <- sum(tt)
		tmps <- matrix(NA, nrow = nn, ncol = n)
		t <- rep(NA, n)
		S0 <- t
		t[!tt]  <- tmp$surv0[, 1]    # times
		S0[!tt] <- tmp$surv0[, 2]    # baseline survival function

		for (i in 1:nn) {
			tmps[i,] <- S0^exp(lp[i])
		}

		xsurv[tt, ] <- tmps  # individual by times matrix of xvalidated survival curves
	}

	res <- survinterp(y, t(xsurv))

# get observed and fitted counts
#	tmp <- eval.surv.fit(y, event, res$xsurv, breaks = breaks)

# return crossvalidated linear predictor and cross validated and
# interpolated survival curves
   k<-order(y)
   C<-matrix(c(C))
   attributes(C)<-list(grp=grp)
	return(list(
		xval.lp = C,
		times=y[k],
		event.sort=event[k],
#		obs.count = tmp$observed,
#		exp.count = tmp$expected,
		xval.surv = res$xsurv
	))
}
