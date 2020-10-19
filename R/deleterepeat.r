DeleteRepeat<-function (X, y, event = NULL, method = HGmultc, rounds = 1,xvalid = TRUE
                        ,fold = NULL, trace = FALSE, weights = rep(1, nrow(X)), ...)
{
    weights <- abs(weights)
    nn <- sum(weights > 0)
    weights <- (weights/sum(weights)) * nn
    type <- "class"
    xfv <- matrix(0, nrow = length(y), ncol = rounds)
    xvsurv<-list()
    k<-2:(nrow(X)+1)
    if(xvalid) grp<-  matrix(0, nrow = length(y), ncol = rounds)  # added
    errors <- rep(-1, rounds)
    when <- rep(-1, ncol(X))
    for (i in 1:rounds) {
    if(!any(when<0)){
    cat("stopping after round",i-1,": no variables left","\n")
    break
    }
        if (trace)
            cat("DeleteRepeat round ", i, "\n")
        if (is.null(event))
             tt <- method(X[, when < 0], y, weights,...)  #added w(=weights) here
        else {tt <- method(X[, when < 0], y, event, ...)
            xvsurv[[i]]<-cbind(tt$surv0[,1],tt$surv0[,3],tt$survall[,-1] )
            colnames(xvsurv[[i]])[1:2]<-c("time","event")
            }
        if (xvalid) {
            fv <- xvalidate(X[, when < 0], y, event, method, fold,
                trace = trace, weights, ...)
            if(is.list(fv)){
            xvsurv[[i]]<-fv$xval.surv
            fv<-fv[[1]]}
            errors[i] <- error(tt, X[, when < 0], y, weights, fv)
        }
        else {
            errors[i] <- error(tt, X[, when < 0], y, weights)
            fv <- predict(tt, X[, when < 0])
        }
        ind <- which(when < 0)[which.genes(tt)]
        when[ind] <- i
        xfv[, i] <- fv
        if(xvalid)grp[,i] <- attr(fv,"grp")
    }
    if(xvalid)
    return(list(genes.when = when, error.rates = errors, xvfv = xfv, grp=grp, xvsurv=xvsurv))
    else
    return(list(genes.when = when, error.rates = errors, xvfv = xfv, xvsurv=xvsurv))
}
