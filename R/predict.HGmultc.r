predict.HGmultc <- function (object, newdata, type = c("class", "prob"), ...)
{
	type <- match.arg(type)

	if (missing(newdata))
		prob <- object$P
	else
		prob <- HGIprob(cbind(1, newdata), object$beta, object$S)

	if (type == "prob")
		return(prob)
	else
		return(apply(prob, 1, which.max))
}

