test.HGordcat <-
  function(testing   = TRUE,
           test.name = "HGordcat"
           )
  {
    set.seed(1)

    ## EXAMPLE - Data not generated by model but still works
    n      <- 40
    x      <- matrix(rnorm(n * n), nrow = n)
    p      <- matrix(0, nrow = n, ncol = 3)
    vec    <- quantile(x[, 1], probs = c(0, 1/3, 2/3, 1))
    vec[1] <- vec[1] - 0.1
    vec[4] <- vec[4] + 0.1
    y      <- as.numeric(cut(x[, 1], vec))
    try    <- HGordcat(x, y)

    ## Check results
    test.compare(try,
                 truth,
                 test.name,
                 testing
                 )
  }
