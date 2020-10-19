"checkArgs" <- function(X, y, weights, lambda, bbess, kbess, adj.prior, Call) {
	Objs <- ls()
	Objs <- Objs[-match("Call", Objs)]
	FunName <- deparse(Call[[1]])
	ArgNames <- names(Call)[-1]
	Objs <- Objs[-(match(ArgNames, Objs))]

#	cat("Call:", as.character(Call), "\n")
#	cat("match.call:", as.character(match.call(call = sys.call(sys.parent()))), "\n")

	flags <- list()
	if (!is.null(ArgNames)) {
		for (i in 2:length(Call)) {
			if (inherits(try(eval(Call[[i]]), TRUE), "try-error")) {
				flags[[ArgNames[i - 1]]] <- FALSE
			} else {
				flags[[ArgNames[i - 1]]] <- TRUE
			}
		}

		if (length(Objs) > 0) {
			for (i in 1:length(Objs)) {
				flags[[Objs[i]]] <- TRUE
			}
		}
	}

	ErrMsg <- c()

	if (missing(X))
		ErrMsg <- c(ErrMsg, "'X' is missing.")
	else if (!is.matrix(X))
		ErrMsg <- c(ErrMsg, "'X' must be a numeric matrix.")
	else if (is.matrix(X) && !is.numeric(X))
		ErrMsg <- c(ErrMsg, "'X' must be a numeric matrix.")

	if (missing(y))
		ErrMsg <- c(ErrMsg, "'y' is missing.")
	else if (!is.vector(y, mode = "numeric"))
		ErrMsg <- c(ErrMsg, "'y' must be a numeric vector.")

	if (!missing(X) && !missing(y)) {
		if (is.matrix(X) && is.vector(y, mode = "numeric")) {
			if (length(y) != nrow(X))
				ErrMsg <- c(ErrMsg, "Number of elements in 'y' must be equal to the number of rows in 'X'.")
			if (any(is.na(y)) || any(is.na(X)))
				ErrMsg <- c(ErrMsg, "No missing values allowed in 'y' or 'X'.")

			if (missing(weights))
				ErrMsg <- c(ErrMsg, "Weights are missing.")
			else {
				if (any(weights <= 0))
					ErrMsg <- c(ErrMsg, "Weights must be positive")
				if (length(weights) != length(y))
					ErrMsg <- c(ErrMsg, "The number of weights is not equal to the number of elements of 'y'.")
				if (length(weights) != nrow(X))
					ErrMsg <- c(ErrMsg, "The number of weights is not equal to the number of rows in 'X'.")
			}
		}
	}

	is.scalar <- function(z) {
		if (is.vector(z, mode = "numeric") && length(z) == 1)
			z <- TRUE
		else z <- FALSE
			z
	}

	if (!is.scalar(bbess))
		ErrMsg <- c(ErrMsg, "'bbess' must be a scalar.")
	else
		if (bbess < 1 | bbess > 1e+07)
			ErrMsg <- c(ErrMsg, "'bbess' is outside the allowable range [1, 1e+07].")

	if (!is.scalar(kbess))
		ErrMsg <- c(ErrMsg, "'kbess' must be a scalar.")
	else
		if (kbess < 0 | kbess > 1)
			ErrMsg <- c(ErrMsg, "'kbess' is outside the allowable range [0, 1].")

	switch (FunName,
		TestFn1 = {
			if (!is.scalar(lambda))
				ErrMsg <- c(ErrMsg, "'lambda' must be a positive scalar.")
			else if (is.scalar(lambda) && lambda <= 0)
				ErrMsg <- c(ErrMsg, "'lambda' must be a positive scalar.")

			if (!missing(y) && is.vector(y, mode = "numeric") && !any(is.na(y))) {
				if (any(y < 1))
					ErrMsg <- c(ErrMsg, "'y' must be a real-valued vector of class labels 1, 2, ...,C.")
				if (length(unique(y)) == 1)
					ErrMsg <- c(ErrMsg, "'y' must be a real-valued vector of class labels 1, 2, ...,C.")
			}

			if (any(is.na(adj.prior)))
				ErrMsg <- c(ErrMsg, "'adj.prior' must be either TRUE or FALSE.")
			else {
				if (is.vector(adj.prior, mode = "logical") && length(adj.prior) != 1)
					ErrMsg <- c(ErrMsg, "'adj.prior' must be either TRUE or FALSE.")
				if (!is.vector(adj.prior, mode = "logical"))
					ErrMsg <- c(ErrMsg, "'adj.prior' must be either TRUE or FALSE.")
			}
		}
	)

	if (length(ErrMsg) == 0)
		return(NULL)
	else
		return(ErrMsg)
}

#Die1 <- function(whinge) {
#	if(length(whinge) == 1)
#	stop(whinge, call. = FALSE)
#	else {
#	for(i in seq(whinge))
#			warning(whinge[i], call. = FALSE)
#	stop("Multiple errors encountered.", call. = FALSE)
#	}
#}

Halt <- function(Msg) {
	lMsg <- length(Msg)
	headline <- paste(lMsg, ngettext(lMsg, "error", "errors"), "occurred:\n")
	AllMsg <- paste(headline, paste(Msg, collapse = "\n"), collapse = "\n")
	stop(AllMsg, call. = FALSE)
}
