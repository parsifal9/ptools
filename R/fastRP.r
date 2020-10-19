fastRP <- function (X, Y, weights = rep(1, nrow(X)), offset = NULL,
	RPmethod = "class", x = FALSE, y = TRUE, parms, control, cost, ...)
{
	if (is.null(colnames(X)))
		colnames(X) <- 1:ncol(X)  # added patch

	call <- match.call()

	if (is.null(weights))
		weights <- rep(1, nrow(X))  # added

	nobs <- nrow(X)
	nvar <- ncol(X)
	method <- RPmethod       #added

	if (missing(method)) {
		if (is.factor(Y) || is.character(Y) || is.integer(Y))
			method <- "class"
		else if (inherits(Y, "Surv"))
			method <- "exp"
		else if (is.matrix(Y))
			method <- "poisson"
		else
			method <- "anova"
	}

	if (is.list(method)) {
		mlist <- method
		method <- "user"

		if (missing(parms))
			init <- mlist$init(Y, offset, wt = weights)
		else
			init <- mlist$init(Y, offset, parms, weights)

		method.int <- 4
		keep <- rpartcallback(mlist, nobs, init)
	}
	else {
		method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))

		if (is.na(method.int))
			stop("Invalid method")

		method <- c("anova", "poisson", "class", "exp")[method.int]

		if (method.int == 4)
			method.int <- 2

		if (!exists("getNamespace")) {
			init.fun <- get(paste("rpart", method, sep = "."))
		}
		else {
			init.fun <- get(paste("rpart", method, sep = "."),
				getNamespace("rpart"))
		}

		if (missing(parms))
			init <- init.fun(Y, offset, , weights)
		else
			init <- init.fun(Y, offset, parms, weights)
	}

	Y <- init$y

	if (is.matrix(X) && is.numeric(X)) {
		xlevels <- NULL
		isord <- rep(FALSE, nvar)
	}
	else {
		X <- data.frame(X)
		isf <- sapply(X, is.factor)
		xlevels <- lapply(X, levels)
		names(xlevels) <- names(X)
		xlevels <- xlevels[isf]
		isord <- sapply(X, is.ordered)
		X <- sapply(X, as.numeric)
	}

	cats <- rep(0, nvar)

	if (!is.null(xlevels))
		cats[match(names(xlevels), dimnames(X)[[2]])] <-
			unlist(lapply(xlevels, length))

	extraArgs <- list(...)

	if (length(extraArgs)) {
		controlargs <- names(formals(rpart.control))
		indx <- match(names(extraArgs), controlargs, nomatch = 0)

		if (any(indx == 0))
			stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
	}

	controls <- rpart.control(...)

	if (!missing(control))
		controls[names(control)] <- control

	xval <- controls$xval

	if (is.null(xval) || (length(xval) == 1 && xval == 0) || method == "user") {
		xgroups <- 0
		xval <- 0
	}
	else if (length(xval) == 1) {
		xgroups <- sample(rep(1:xval, length = nobs), nobs, replace = FALSE)
	}
	else if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
	}
	else
		stop("Wrong length for xval")

	if (missing(cost))
		cost <- rep(1, nvar)
	else {
		if (length(cost) != nvar)
			stop("Cost vector is the wrong length")

		if (any(cost <= 0))
			stop("Cost vector must be positive")
	}

	rpfit <- .C("s_to_rp",
		n = as.integer(nobs),
		nvarx = as.integer(nvar),
		ncat = as.integer(cats * (!isord)),
		method = as.integer(method.int),
		as.double(unlist(controls)),
		parms = as.double(unlist(init$parms)),
		as.integer(xval),
		as.integer(xgroups),
		as.double(t(init$y)),
		as.double(X),
		as.integer(!is.finite(X)),
		error = character(1),
		wt = as.double(weights),
		as.integer(init$numy),
		as.double(cost),
		NAOK = TRUE, PACKAGE = "rpart")

	if (rpfit$n == -1)
		stop(rpfit$error)

	nodes <- rpfit$n
	nsplit <- rpfit$nvarx
	numcp <- rpfit$method
	ncat <- rpfit$ncat[1]
	numresp <- init$numresp

	if (nsplit == 0)
		xval <- 0

	cpcol <- if (xval > 0 && nsplit > 0) 5 else 3

	if (ncat == 0)
		catmat <- 0
	else
		catmat <- matrix(integer(1), ncat, max(cats))

	rp <- .C("s_to_rp2",
		as.integer(nobs),
		as.integer(nsplit),
		as.integer(nodes),
		as.integer(ncat),
		as.integer(cats * (!isord)),
		as.integer(max(cats)),
		as.integer(xval),
		which = integer(nobs),
		cptable = matrix(double(numcp * cpcol), nrow = cpcol),
		dsplit = matrix(double(1), nsplit, 3),
		isplit = matrix(integer(1), nsplit, 3),
		csplit = catmat,
		dnode = matrix(double(1), nodes, 3 + numresp),
		inode = matrix(integer(1), nodes, 6),
		PACKAGE = "rpart")

	tname <- c("<leaf>", dimnames(X)[[2]])

	if (cpcol == 3)
		temp <- c("CP", "nsplit", "rel error")
	else
		temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")

	dimnames(rp$cptable) <- list(temp, 1:numcp)

	dn1 <- if (nsplit == 0)
		character(0)
	else
		tname[rp$isplit[, 1] + 1]

	splits <- matrix(c(rp$isplit[, 2:3], rp$dsplit), ncol = 5,
		dimnames = list(dn1, c("count", "ncat", "improve", "index", "adj")))

	index <- rp$inode[, 2]
	nadd <- sum(isord[rp$isplit[, 1]])

	if (nadd > 0) {
		newc <- matrix(integer(1), nadd, max(cats))
		cvar <- rp$isplit[, 1]
		indx <- isord[cvar]
		cdir <- splits[indx, 2]
		ccut <- floor(splits[indx, 4])
		splits[indx, 2] <- cats[cvar[indx]]
		splits[indx, 4] <- ncat + 1:nadd

		for (i in 1:nadd) {
			newc[i, 1:(cats[(cvar[indx])[i]])] <- -1 * as.integer(cdir[i])
			newc[i, 1:ccut[i]] <- as.integer(cdir[i])
		}

		if (ncat == 0)
			catmat <- newc
		else catmat <- rbind(rp$csplit, newc)
			ncat <- ncat + nadd
	}
	else
		catmat <- rp$csplit

	if (nsplit == 0) {
		frame <- data.frame(
			row.names = 1,
			var = "<leaf>",
			n = rp$inode[, 5],
			wt = rp$dnode[, 3],
			dev = rp$dnode[, 1],
			yval = rp$dnode[, 4],
			complexity = rp$dnode[, 2],
			ncompete = pmax(0, rp$inode[, 3] - 1),
			nsurrogate = rp$inode[, 4]
		)
	}
	else {
		temp <- ifelse(index == 0, 1, index)
		svar <- ifelse(index == 0, 0, rp$isplit[temp, 1])

		frame <- data.frame(
			row.names = rp$inode[, 1],
			var = factor(svar, 0:ncol(X), tname),
			n = rp$inode[, 5],
			wt = rp$dnode[, 3],
			dev = rp$dnode[, 1],
			yval = rp$dnode[, 4],
			complexity = rp$dnode[, 2],
			ncompete = pmax(0, rp$inode[, 3] - 1),
			nsurrogate = rp$inode[, 4]
		)
	}

	if (method.int == 3) {
		numclass <- init$numresp - 1
		temp <- rp$dnode[, -(1:4)] %*%
			diag(init$parms$prior * sum(init$counts)/pmax(1, init$counts))

		yprob <- temp / rowSums(temp)
		yval2 <- matrix(rp$dnode[, -(1:3)], ncol = numclass + 1)

		frame$yval2 <- cbind(yval2, yprob)
	}
	else if (init$numresp > 1)
		frame$yval2 <- rp$dnode[, -(1:3)]

	if (is.null(init$summary))
		stop("Initialization routine is missing the summary function")

	if (is.null(init$print))
		functions <- list(summary = init$summary)
	else
		functions <- list(summary = init$summary, print = init$print)

	if (!is.null(init$text))
		functions <- c(functions, list(text = init$text))

	if (method == "user")
		functions <- c(functions, mlist)

	where <- rp$which
	names(where) <- row.names(X)

	if (nsplit == 0) {
		ans <- list(frame = frame, where = where, call = call,
			cptable = t(rp$cptable), method = method, parms = init$parms,
			control = controls, functions = functions)
	}
	else {
		ans <- list(frame = frame, where = where, call = call,
			cptable = t(rp$cptable), splits = splits, method = method,
			parms = init$parms, control = controls, functions = functions)
	}

	if (ncat > 0)
		ans$csplit <- catmat + 2

	if (y)
		ans$y <- Y

	if (x) {
		ans$x <- X
		ans$wt <- weights
	}

	ans$ordered <- isord

	if (!is.null(xlevels))
		attr(ans, "xlevels") <- xlevels

	if (method == "class")
		attr(ans, "ylevels") <- init$ylevels

	ans$gene.names <- colnames(X)
	class(ans) <- c("fastRP", "rpart")

	ans
}

