substructure <- function(X, y, maxrounds = 5) {
  dl <- dlda(X, y)
  m1 <- sdda(dl, X, y)
  nn <- which.genes(m1)

  p <- ncol(X)

  k <- length(nn)
  genes <- vector("list", k)

  for (i in 1:k) {
    st <- rep(FALSE, p)
    st[nn[-i]] <- TRUE

    nv <- rep(FALSE, p)
    nv[nn[i]] <- TRUE

    grps <- vector("list", maxrounds + 1)
    grps[[1]] <- nn[i]

    for (j in 1:maxrounds) {
      mm <- sdda(dl, X, y, start = st, never = nv)
      tt <- mm$S & !st

      if (any(tt)) {
        grps[[j + 1]] <- which(tt)
        nv[tt] <- TRUE
      }
      else {
      	break
      }
    }
    genes[[i]] <- grps
  }

  class(genes) <- "substructure"
  return(genes)
}

toplevel <- function(obj) sapply(obj, function(x) x[[1]])

genegroup <- function(obj,i)  unlist(obj[[i]])

print.substructure <- function(x, gene.names = NULL, file = "",
	indent = "    ", ...)
{
	if (is.null(gene.names))
		gene.names <- 1:max(unlist(x))

	if (file != "") cat("Print of substructure object",
		deparse(substitute(x)), "\n", file = file)

	OL <- function(...) cat(..., "\n", file = file, append = TRUE)
	OL("First round model has genes", gene.names[toplevel(x)])

	for (i in 1:length(x)) {
		tt <- x[[i]]
		OL("Groups associated with", gene.names[tt[[1]]])

		for (j in 2:length(tt))
			OL(indent, gene.names[tt[[j]]])
	}

	invisible(x)
}

substructurePlot <- function(X, y, obj, i = 0, lab = yclass, col = yclass + 1,
	...)
{
	if (i == 0) {
		PairsPlot(X[, toplevel(obj)], lab = lab, col = col, ...)
	} else {
		if (all(sapply(obj[[i]], length) <= 1)) {
			PairsPlot(X[, genegroup(obj, i)], lab = lab, col = col, ...)
		}
		else {
#			require(mva)
			tt <- obj[[i]]
			Y <- X[, tt[[1]]]
			nms <- colnames(X)[tt[[1]]]

			for (j in 2:length(tt)) {
				if (length(tt[[j]]) > 1) {
					cxy <- cancor(X[, tt[[1]]], X[, unlist(tt[[j]])])
					Y <- cbind(Y, X[, unlist(tt[[j]])] %*% cxy$ycoef[, 1])
					nms <- c(nms, paste(colnames(X)[unlist(tt[[j]])], sep = "",
						collapse = "."))
				} else if (length(tt[[j]]) == 1) {
					Y <- cbind(Y, X[, tt[[j]]])
					nms <- c(nms, colnames(X)[tt[[j]]])
				}
			}

			colnames(Y) <- nms
			PairsPlot(Y, lab = lab, col = col, ...)
		}
	}

	invisible()
}
