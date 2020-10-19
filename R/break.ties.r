# break ties function - use before running HGcoxreg
break.ties <- function(x, sd.divisor = 10)
{
	i <- order(x)
	d <- abs(diff(x[i]))
	sd <- sqrt(min(d[d != 0]))

	if (any(d == 0)) {
		index <- (2:length(x))[d == 0]
		x[i][index] <- x[i][index] + rnorm(length(index), sd = sd / sd.divisor)
	}

	x
}
