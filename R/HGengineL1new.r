HGengineL1new<-function (X, y, event, weights = rep(1, nrow(X)), bbess = 1e+07,
    kbess = 0, scale, bhat, Le, DLde, D2Lde, HGsc, sparsity.prior = "NG",
    no.prior = NULL, Offset = NULL, Offset.par.update = NULL,
    a = NULL, control = HGcontrol())
{
    eps<-1e-6   # Lange and sinsheimer fix for weights
    zpar<-0.5/.Machine$double.eps
    scale.fix <- FALSE
    if (scale < 0) {
        scale.fix <- TRUE
        scale <- abs(scale)
    }
    tolc <- control$tolc
    HGcexp <- match.prior(sparsity.prior)
    n <- nrow(X)
    p <- ncol(X)
    weights <- abs(weights)
    nn <- sum(weights > 0)
    weights <- (weights/sum(weights)) * nn
    offset <- NULL
    if (!is.null(Offset))
        offset <- Offset(y, event, a)
    S <- bhat != 0
    S[no.prior] <- TRUE
    R <- HGeta(X, bhat, S, offset, eta.lim = control$eta.lim)
    #weights<-pmin(1/abs(y-R),zpar)  # added for L1 regression
    weights<-1/(eps+abs(y-R))
    Q1 <- HGpost(y, event, bhat, bhat, S, R, weights, bbess,
        kbess, Le, scale, HGcexp = HGcexp, no.prior)
    beta <- bhat
    j <- 1
    repeat {
        Q0 <- Q1
        repeat {
            delta <- HGdirection(X, y, event, beta, bhat, S,
                R, weights, bbess, kbess, DLde, D2Lde, scale,
                HGcexp = HGcexp, no.prior)
            alpha <- 1
            nc <- 0
            repeat {
                btmp <- beta
                tt <- S
                btmp[tt] <- btmp[tt] + alpha * delta
                R <- HGeta(X, btmp, S, offset,eta.lim = control$eta.lim)
                Q2 <- HGpost(y, event, btmp, bhat, S, R, weights,
                  bbess, kbess, Le, scale, HGcexp = HGcexp, no.prior)
                if (Q2 > Q1 - tolc)
                  break
                nc <- nc + 1
                if (nc > 100)
                  stop("Error in line search. No upward step could be found")
                alpha <- alpha/2
            }
            beta <- btmp
            if (abs(Q2 - Q1) < tolc)
                break
            Q1 <- Q2
            if (alpha * max(abs(delta)) < control$tolerance) {
                break
            }
        }
        j <- j + 1
        if (exists("Qlog"))
            Qlog <<- c(Qlog, Q1)
        S <- (abs(beta) > control$epsilon)
        S[no.prior] <- TRUE
        bhat <- beta
        R<-HGeta(X,bhat,S,offset,eta.lim = control$eta.lim)     # added for L1 regression
       # weights<-pmin(1/abs(y-R),zpar)  # added for L1 regression
        weights<-1/(eps+abs(y-R))
        #cat("weights","\n")
        #cat(weights,"\n")
        #lval<- -sum(abs(y-R))-(2/bbess)^0.5*sum(abs(bhat))
        #cat("criterion",lval,"\n")
        if (abs(Q1 - Q0) < control$tolerance)
            break
        if (!scale.fix)
            scale <- HGsc(R, y, event, weights, scale)
        if (!is.null(Offset)) {
            a <- Offset.par.update(y, event, a, R, weights)
            offset <- Offset(y, event, a = a)
        }
        R <- HGeta(X, bhat, S, offset,eta.lim = control$eta.lim)
        Q1 <- HGpost(y, event, bhat, bhat, S, R, weights, bbess,
            kbess, Le, scale, HGcexp = HGcexp, no.prior)
        if (j > control$maxit)
            break
    }
    zz <- list(beta = bhat, S = S, eta = R, scale = scale, a = a)
    class(zz) <- c("HGengine", "HG")
    zz
}


