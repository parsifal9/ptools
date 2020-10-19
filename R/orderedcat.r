orderedcat <- function (x, y, G, weights)
{
    n <- nrow(x)
    G1 <- G - 1
    tmp <- matrix(rep(weights, G1), nrow = n)
    weights <- as.vector(t(tmp))
    xs <- matrix(0, nrow = n * G1, ncol = G1 + ncol(x))
    I <- diag(G1)
    i <- rep(1:n, rep(G1, n))
    xs[, G:ncol(xs)] <- x[i, ]

    for (j in 1:G1) {
        xs[, j] <- rep(I[, j], n)
    }

    I <- matrix(0, nrow = n, ncol = G)

    for (j in 1:n) {
        I[j, y[j]] <- 1
    }

    I <- apply(I, 1, cumsum)
    ys <- as.vector(I[1:G1, ])

    list(xs = xs, ys = ys, weights = weights)
}

