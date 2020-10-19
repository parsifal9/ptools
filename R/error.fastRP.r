error.fastRP <- function (object, X, y, weights, fv, ...)
{
	if (object$method == "class") {
		y <- as.numeric(y)

		if (missing(fv)) {
			b <- tapply(weights, list(y, predict.fastRP(object, X, type = "class")),
				simplify = TRUE, FUN = "sum")
			e <- 1 - sum(diag(b)) / sum(weights)
		}
		else {
			e <- sum(weights * (y != fv)) / sum(weights)
		}
	}
	else {
		if (missing(fv)) {
		#	browser()
			fv <- predict.fastRP(object, X, type = "anova")
			r <- (y - fv)
		}
		else {
			r <- (y - fv)
			e <- mean(r * r * weights)
		}
	}

	e
}
