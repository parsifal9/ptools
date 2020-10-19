"which.genes.HGsurv" <-
function (object, ...)
{
    S <- object$S
    if (!is.null(dim(S)))
        S <- apply(S, 1, any)
    return(which(S)-1)
}
