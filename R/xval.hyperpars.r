xval.hyperpars <- function(x, y, method = HGmultc, fold = NULL, foldbk = fold,
	trace = TRUE, etol=.05, bbess = 1e7, kbess = seq(0, 1, .1),
	b0sc = 5, ...)
{
# cross validates choice of hyperparameters in HG functions
	ind <- any(names(formals(method)) == "b0sc")

	n <- nrow(x)
	xvfv <- rep(NA,n)
	if (is.null(fold)) {
		grp <- 1:n
	}
	else {
		grp <- rep(1:fold, n/fold + 1)[1:n]
		grp[order(as.numeric(y) + runif(n))] <- grp
	}

	if (ind) {
		details <- matrix(0, nrow = fold, ncol = 5)
		colnames(details) <- c("fold", "kbess", "bbess", "b0sc", "xval-error")
		rownames(details) <- rep("", fold)

		for (i in sort(unique(grp))) {
			tt <- grp == i
			cat("Processing fold", i, "\n")

			res <- optim.hyperpars(x[!tt, ], y[!tt], method = method,
				fold = foldbk, bbess = bbess,
				kbess = kbess, b0sc = b0sc, xvalid = TRUE, trace = trace, ...)

			# get largest b, and smallest k and b0sc within etol of minimum xvalid error
			a <- res$xvtable
			ii <- which(a <= (1+etol)*min(a))
			da <- dim(a)
			nda <- prod(da)
			i1 <- as.numeric(gl(da[1], 1, nda))
			i2 <- as.numeric(gl(da[2], da[1], nda))
			i3 <- as.numeric(gl(da[3], da[1] * da[2], nda))

		#	b <- min(bbess[i2[ii]])
		#	k <- min(kbess[i1[ii]])
		#	b0 <- min(b0sc[i3[ii]])
      j<-order(i1[ii],i3[ii],-i2[ii])
      l<-j[1]
      k<-kbess[i1[ii][l]]
      b0<-b0sc[i3[ii][l]]
      b<-bbess[i2[ii][l]]
      amin<-a[ii][l]
      
			details[i, ] <- c(i, k, b, b0, amin)
			res <- method(x[!tt, ], y[!tt], bbess = b, kbess = k, b0sc = b0, ...)
			xvfv[tt] <- predict(res, x[tt, , drop = FALSE], type = "class")
		}
	}
	else {
		details <- matrix(0, nrow = fold, ncol = 4)
		colnames(details) <- c("fold", "kbess", "bbess", "xval-error")
		rownames(details) <- rep("", fold)

		for (i in sort(unique(grp))) {
			tt <- grp == i
			cat("Processing fold", i, "\n")

			res <- optim.hyperpars(x[!tt, ], y[!tt], method = method,
				fold = foldbk, bbess = bbess,
				kbess = kbess, xvalid = TRUE, trace = trace,...)

			# get largest b and smallest k within prop of minimum xvalid error
			a <- res$xvtable
			da <- dim(a)
			ii <- which(a <= (1+etol)*min(a))
			i1 <- as.numeric(gl(da[1], 1, da[2] * da[1]))
			i2 <- as.numeric(gl(da[2], da[1]))

		#	b <- max(bbess[i2[ii]])
		#	k <- min(kbess[i1[ii]])
      j<-order(i1[ii],-i2[ii])
      l<-j[1]
      k<-kbess[i1[ii][l]]
      b<-bbess[i2[ii][l]]
      amin<-a[ii][l]
			details[i, ] <- c(i, k, b, amin)
			res <- method(x[!tt, ], y[!tt], bbess = b, kbess = k, ...)
			xvfv[tt] <- predict(res, x[tt, , drop = FALSE], type = "class")
		}
	}

	list(xvfv = xvfv, grp = grp, details = details)
}
