xvalidate<-function (X, y, event = NULL, method, fold = NULL, trace = FALSE,
    weights = rep(1, nrow(X)), ...)
{
    if (is.null(colnames(X)))
        colnames(X) <- 1:ncol(X)
    type <- "class"
    n <- nrow(X)
    C <- rep(0, n)
    if (!is.null(event)) {
        res <- xval.survival(X, y, event, method = method, fold = fold,
            trace = trace, weights = weights,...)
        return(list(xval.lp = res[[1]], xval.surv = cbind(res[[2]],res[[3]],res[[4]])))
    }
    else {
        if (is.null(fold)) {
            grp <- 1:n
        }
        else {
            grp <- rep(1:fold, n/fold + 1)[1:n]
            grp[order(as.numeric(y) + runif(n))] <- grp
        }
        for (i in sort(unique(grp))) {
            if (trace)
                cat("X-validate group", i, "of", max(grp), "\n")
            tt <- (grp == i)
            if (is.null(event))
                tmp <- method(X[!tt, , drop = FALSE], y[!tt],
                  weights = weights[!tt], ...)
            else tmp <- method(X[!tt, , drop = FALSE], y[!tt],
                event[!tt], weights = weights[!tt], ...)
            if (class(tmp)[1] == "fastRP") {
                type <- "class"
                if (tmp$method == "anova")
                  type <- "vector"
            }
            C[tt] <- predict(tmp, X[tt, , drop = FALSE], type = type)
        }
        if (is.null(type)) {
            names(C) <- make.names(grp)
            return(C)
        }
        else {
        C<-matrix(C,nrow=n,byrow=T)
        attributes(C)<-list(grp=grp) # attributes added here
        return(C)
        }
    }
}
