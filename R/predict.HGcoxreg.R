"predict.HGcoxreg" <- function(object, newdata, type = type, ...)
{
# at the moment simply return linear predictor
	if (length(object$beta) == 1)
		object$coeffs <- as.matrix(object$beta)

#	coeffs<-object$beta[object$S]
	fv <- as.matrix(newdata) %*% object$beta
	fv
}
