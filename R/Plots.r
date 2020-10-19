PairsPlot <- function(x, y, lab = 1, col = 1) {
	if (missing(y)) {
		y <- x
		yflg <- FALSE
	} else {
		yflg <-  TRUE
	}

	n <- ncol(x)
	m <- ncol(y)
	gap <- 0.2
	botleft <- 5.1
	topright <- 2.1
	op <- par(mfrow = c(n, m), mar = rep(gap, 4),
		oma = c(botleft, botleft, topright, topright))

	for (i in 1:n) {
		for (j in 1:m) {
			plot(y[, j], x[, i], type = "n", axes = FALSE)
			box()

			if (yflg || i != j)
				points(y[, j], x[, i], pch = lab, col = col)

			if (i == n) {
				axis(1, xpd = FALSE)
				mtext(colnames(y)[j], side = 1, line = 3.5)
			}

			if (j == 1) {
				axis(2, xpd = FALSE)
				mtext(colnames(x)[i], side = 2, line = 3.5)
			}
		}
	}

	par(op)
}
