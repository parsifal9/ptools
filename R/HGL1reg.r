 # version using Hunter and Lange MM approximation to abs
 # source the files in Rchip3.0.2 and rename this to HGL1reg to use error and
 # predict functions in RCHIP
HGL1reg<-function (x, f, weights = rep(1, nrow(x)), sparsity.prior = "NG",
    bbess = 1e7, kbess = 1, b0sc = 5, scale = -1, initb = NULL,
    no.prior = NULL)
{
    if(length(initb)==ncol(x))initb<-c(initb,rnorm(1))
              lambda <- 1
    if (nrow(x) > ncol(x)) {
        lambda <- 1e-04
        b0sc <- -1
    }
    x <- cbind(rep(1, nrow(x)), x)
    if (is.null(initb)) {
        bhat <- HGinit.g(x, f, eps = 0, lb = lambda, weights,
            b0sc)
        nl <- max(abs(bhat[-1]))
        if (b0sc < 0)
            b0sc <- nl
        ysc <- b0sc/nl
        bhat <- bhat * ysc
    }
    else {
        bhat <- initb
        ysc <- 1
        }

    res <- HGengineL1new(x, f * ysc, event = NULL, weights = weights,
        bbess = bbess, kbess = kbess, scale = scale, bhat = bhat,
        Le = Le.g, DLde = DLde.g, D2Lde = D2Lde.g, HGsc = HGsc.g,
        sparsity.prior = sparsity.prior, no.prior = no.prior,
        control = HGcontrol(tolc = 1e-05,epsilon=1e-4,tolerance=1e-8,maxit=300))
    sigma2 <- res$scale/ysc^2
    if (scale < 0)
        sigma2 <- abs(scale)
    zz <- list(beta = res$beta/ysc, S = res$S, fv = res$eta/ysc,
        varids = which(res$S) - 1, sigma2 = sigma2)
    class(zz) <- c("HGL1reg", "HG")
    zz
}

