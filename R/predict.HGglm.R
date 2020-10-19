"oldpred" <- function(object, newdata, type = type, ...)
{
	if (length(object$coeffs) == 1)
		object$coeffs <- as.matrix(object$coeffs)

	coeffs <- object$beta[object$S]

	fv <- fvalfn(cbind(rep(1, nrow(newdata)), newdata), object$S, coeffs)
	fv
}

"predict.HGglm" <- function(object, newdata, type = type, ...)
{
	if (length(object$coeffs) == 1)
		object$coeffs <- as.matrix(object$coeffs)

	coeffs <- object$beta[object$S]
	n <- nrow(newdata)
	lp <- (cbind(rep(1, n), newdata)[, object$S]) %*% coeffs

	model <- object$model
	if (model == "N")   fv <- lp
	if (model == "B")   fv <- exp(lp)/(1+exp(lp))
	if (model == "P")   fv <- exp(lp)
	if (model == "G")   fv <- 1/lp
	if (model == "IG")  fv <- 1/lp^2
	if (model == "Own") fv <- fvalfn(lp)

	fv
}
