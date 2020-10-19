eval.surv.fit <- function(times, event, surv,
	breaks = quantile(times, probs = c(.25, .5, .75)))
{
# times is vector of new survival times n by 1
# event is new censoring vector n by 1
# ts is suppport of survival functions ns by 1
# surv is an ns by n matrix of survival curves
# over ts relevant for times and event
# works out nclass by nclass table of observed and expected deaths
# the 1/nclass quantiles of survival times are used to define the classes
# breaks will be used in cut by adding 0 and max(times) to get time bins
# could use quantiles of times eg breaks=quantiles(times,probs=c(1/3,2/3))

	td <- times
	ts <- times
	S <- surv
	#probs <- seq(1/nclass,1-1/nclass,1/nclass)
	#tc <- quantile(td,probs=probs)
	#tc <- c(5,10,15,20)
	tc <- breaks
	breaks <- c(0, tc, max(times))
	a <- NULL
	tf <- cut(times,breaks = breaks, labels = FALSE)
	pstar <- surv.prob(ts, S, times)

	for (i in tc) {
		p <- surv.prob(ts, S, i)
		a <- cbind(a, p)
	}

	get.probs <- function(p)
	{
		q <- length(p)
		c(1 - p[1], -diff(p), p[q])
	}

	res <- t(apply(a, 1, get.probs))

	# get expected true values for censored observations
	q <- matrix(0, nrow = 1, ncol = ncol(res))
	ic <- event == 0

	if (any(ic)) {
		q <- res[ic,]
		tfc <- tf[ic]
		pstar1 <- pstar[ic]
		n <- nrow(q)
		qq <- ncol(q)

		for (i in 1:n) {
			k <- tfc[i]

			if (k > 1)
				q[i, 1:(k-1)] <- 0

			if (k < qq) {
				q[i, k] <- (pstar1[i] - sum(q[i, (k+1):qq])) / pstar1[i]
				q[i, (k+1):qq] <- q[i, (k+1):qq] / pstar1[i]
			}
			else {
				q[i,k] <- 1
			}
		}
	}

	#browser()
	#compute true counts in bins
	a <- rep(0, ncol(q))

	if (any(!ic)) {
		aa <- table(tf[!ic])
		ka <- as.numeric(attributes(aa)$dimnames[[1]])
		a[ka] <- as.vector(aa)
	}

	# add contribution from censored observations
	a <- a+colSums(q)

	# get expected numbers in each bin from survival curves
	b <- colSums(res)
	attributes(b) <- NULL
	list(observed = a, expected = b, cut.times = breaks)
}
