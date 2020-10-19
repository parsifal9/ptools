predict.fastRP <- function (object, newdata = list(),
	type = c("class", "prob", "matrix", "vector", "anova"), ...)
{
	if (is.null(colnames(newdata)))
		colnames(newdata) <- 1:ncol(newdata)  # added patch to make it work

    if (!inherits(object, "fastRP"))
        stop("Not legitimate tree")

    mtype <- missing(type)
    type <- match.arg(type)

    if (type == "anova")
    	type <- "vector"  # patch to make predict method more intuitive

    if (!exists("getNamespace")) {
        FUN <- get("pred.rpart")
    }
    else {
        FUN <- get("pred.rpart", getNamespace("rpart"))
    }
    if (missing(newdata)) {
        where <- object$where
    }
    else {
        if (!(is.matrix(newdata) && is.numeric(newdata))) {
            newdata <- data.frame(newdata)
            isf <- sapply(newdata, is.factor)
            xlev <- attr(object, "xlevels")
            for (i in names(newdata)[isf]) {
                newdata[i] <- factor(as.character(newdata[i]),
                  levels = xlev[i])
                newdata[i] <- as.numeric(newdata[i])
            }
            newdata <- as.matrix(newdata)
        }
          #browser()
        where <- FUN(object, newdata)
    }
    frame <- object$frame
    method <- object$method
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if (mtype && nclass > 0)
        type <- "class"
    if (type == "vector" || (type == "matrix" && is.null(frame$yval2))) {
        pred <- frame$yval[where]
        names(pred) <- names(where)
    }
    else if (type == "matrix") {
        pred <- frame$yval2[where, ]
        dimnames(pred) <- list(names(where), NULL)
    }
    else if (type == "class" && nclass > 0) {
        pred <- factor(ylevels[frame$yval[where]], levels = ylevels)
        names(pred) <- names(where)
    }
    else if (type == "prob" && nclass > 0) {
        pred <- frame$yval2[where, 1 + nclass + 1:nclass]
        dimnames(pred) <- list(names(where), ylevels)
    }
    else stop("Invalid prediction for rpart object")
    if (missing(newdata) && !is.null(object$na.action))
        pred <- naresid(object$na.action, pred)
    pred
}

