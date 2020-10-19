which.genes.fastRP <- function (object, ...)
{
    nms <- object$frame[, "var"]
    nms <- unique(nms[nms != "<leaf>"])
    wh <- match(nms, object$genenames)
    #names(wh) <- nms
    #return(wh)
    return(as.numeric(nms) - 1)
}

