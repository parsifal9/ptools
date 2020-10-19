"which.genes.HGmultc" <- function (object, ...) {
    S <- object$S
    if (!is.null(dim(S)))
        S <- apply(S, 1, any)
    S <- S[-1]
    return(which(S))
}
