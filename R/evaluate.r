error.rate <- function(object, X, y, ...) {
  tt <- sum(diag(table(y, predict(object, X, type = "class"))))
  1.0 - as.double(tt) / length(y)
}


permute <- function(X, y, method, B = 100, trace = FALSE, ...) {
  rates <- rep(NA, B)
  tmp <- method(X, y, ...)
  rates[1] <- error.rate(tmp, X, y)

  for (i in 2:B) {
    if (trace) cat("Permutation", i, "of", B, "\n")

    ys <- sample(y)
    tmp <- method(X, ys, ...)
    rates[i] <- error.rate(tmp, X, ys)
  }

  rates
}



