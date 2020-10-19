test.HG <-
  function(fun.names =
           c("HGcoxreg",
             "HGcoxregI",
             "HGgaussian",
             "HGglm",
             "HGmultc",
             "HGordcat",
             "HGsurv"
             )
           )
  {
    runtest <-
      function(fun.name)
        {
          cat("Testing '", fun.name, "' ... ", sep = "")

          test.fun <- get(paste("test", fun.name, sep = "."))

          val <- test.fun()

          ok <- identical(val, TRUE)
          if (ok)
            cat("OK\n")
          else
            cat("FAIL!\n")

          val
        }

    vals <- lapply(as.list(fun.names), runtest)
    names(vals) <- fun.names

    is.true <- function(x) identical(x, TRUE)
    ok <- sapply(vals, is.true)

    if (all(ok))
      {
        cat("All tests passed.\n")
        return(invisible(TRUE))
      }

    test.HG.diffs <<- vals[!ok]
    cat("Not all tests passed. Diagnostic information written to object 'test.HG.diffs'.\n")
    invisible(FALSE)
  }
