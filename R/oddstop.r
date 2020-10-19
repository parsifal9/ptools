oddstop <- function (p, G)
{
    col <- 1
    G <- G - 1
    n <- nrow(p)
    m <- n/G
    pc <- matrix(0, nrow = m, ncol = (G + 1))
    a <- matrix(0, nrow = m, ncol = G + 1)

    for (i in 1:G) {
        a[, i + 1] <- p[seq(i, n, G), col]
    }

    pc[, G + 1] <- a[, G + 1]

    if (G >= 2) {
        for (k in (G + 1):2) {
            pc[, k - 1] <- (a[, k - 1]/a[, k]) * pc[, k] * (1 -
                a[, k])
        }
    }

    pc[, 1] <- 1 - pc %*% rep(1, (G + 1))
    lab <- apply(pc, 1, which.max)

    list(p = pc, lab = lab)
}



