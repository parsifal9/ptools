optim.hyperpars <- function(X, y, event = NULL, method, xvalid = TRUE,
	fold = NULL, trace = FALSE, weights = rep(1, nrow(X)),
	bbess = 1e7, kbess = 0, b0sc = 5, ...)
{
	ind <- any(names(formals(method)) == "b0sc")

	p <- length(b0sc)
	q <- length(bbess)
	r <- length(kbess)

	aa <- rep(NA, p * q * r)

	ii <- 0
	if (ind) {
		b <- list(kbess = kbess, bbess = bbess, b0sc = b0sc)

		for (i in b0sc) {
			for (j in bbess) {
				for (k in kbess) {
					if (trace)
						cat("Computing cross validation error for kbbess =", k,
							"bbess =", j, "b0sc =", i, "\n")

					res <- DeleteRepeat(X, y, event, method = method, rounds = 1,
						xvalid = xvalid, fold = fold, trace = FALSE,
						weights = weights, b0sc = i, bbess = j, kbess = k, ...)

					ii <- ii + 1
					aa[ii] <- res[[2]]
					#cat(i,j,k,aa[ii],"\n")
				}
			}
		}

		a <- array(aa, dim = c(r, q, p), dimnames = b)
	}
	else {
		b <- list(kbess = kbess, bbess = bbess)

		for (j in bbess) {
			for (k in kbess) {
				if (trace)
					cat("Computing cross validation error for kbbess =", k,
						"bbess =", j, "\n")

				res <- DeleteRepeat(X, y, event, method = method, rounds = 1,
					xvalid = xvalid, fold = fold, trace = FALSE, weights = weights,
					bbess = j, kbess = k, ...)

				ii <- ii + 1
				aa[ii] <- res[[2]]
				#cat(j,k,aa[ii],"\n")
			}
		}

		a <- array(aa, dim = c(r, q), dimnames = b)
	}

	list(xvtable = a)
}
