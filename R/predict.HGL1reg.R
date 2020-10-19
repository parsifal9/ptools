"predict.HGL1reg" <- function(object, newdata, type = type, ...)
{
	if (length(object$coeffs) == 1)
		object$coeffs <- as.matrix(object$coeffs)

	fv <- cbind(rep(1, nrow(newdata)), newdata) %*% object$beta
	fv
}
