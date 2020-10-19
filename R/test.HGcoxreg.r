test.HGcoxreg <-
  function(testing   = TRUE,
           test.name = "HGcoxreg"
           )
  {
    set.seed(1)

    n   <- 50
    x   <- matrix(rnorm(n * 200), nrow = n)   # Generate data
    lp  <- 3 + x[, 1]
    a   <- 1                              # Weibull shape parameter 1 - exponential distribution
    tt  <- rweibull(n, a, exp(-lp/a))     # Simulate survival times
    C   <- rep(0:1, length.out = n)       # Alternate observations censored

    set.seed(1)
    try <- HGcoxreg(x, tt, C)             # Fit model

    test.compare(try,
                 truth,
                 test.name,
                 testing
                 )
  }
