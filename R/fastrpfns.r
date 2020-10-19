
serule <- function(tree, ...) {
  UseMethod("serule")
}

serule.rpart <- function(tree, se = 1.0, ...) {
  p <- tree$cptable
  cp0 <- p[,1]
  cp <- sqrt(cp0 * c(Inf, cp0[-length(cp0)]))
  xerror <- p[,4]
  xstd <- p[,5]

  minpos <- min(seq(along = xerror)[xerror == min(xerror)])
  return(cp[xerror < xerror[minpos] + se*xstd[minpos]][1])
}

post.fastRP <- function(tree, title., ...) {
  if (!inherits(tree, "rpart")) stop("Not an rpart object")
  margin <- 0.1
  branch <- 0.2
  all <- TRUE
  digits <- getOption("digits") - 3
  use.n <- FALSE

  dev <- dev.cur()
  if (dev == 1) dev <- 2
  assign(paste(".rpart.parms", dev, sep = "."),
         list(uniform = TRUE, branch = branch, compress = TRUE,
              nspace = branch, minbranch = 0.3),
         envir = .GlobalEnv)
  if (!exists("getNamespace")) {
    FUN <- get("rpartco")
  } else {
    FUN <- get("rpartco",getNamespace("rpart"))
  }
  xy <- FUN(tree)
  xx <- xy$x
  yy <- xy$y
  temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
  temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
  plot(temp1, temp2,
       type = "n", axes = FALSE, xlab = "", ylab = "",
       ...)

  frame <- tree$frame
  node <- as.numeric(row.names(frame))

  if (!exists("getNamespace")) {
    FUN2 <- get("rpart.branch")
  } else {
    FUN2 <- get("rpart.branch",getNamespace("rpart"))
  }

  temp <- FUN2(xx, yy, node, branch)

  if (branch > 0)
    text(xx[1], yy[1], "|")
  lines(c(temp$x), c(temp$y))

  col <- names(frame)
  method <- tree$method
  ylevels <- attr(tree, "ylevels")
  if (!is.null(ylevels)) col <- c(col, ylevels)

  cxy <- c(strwidth("M"),2*strheight("M"))
  if (!is.null(srt <- list(...)$srt) && srt == 90)     cxy <- rev(cxy)
  is.left <- (node%%2 == 0)
  node.left <- node[is.left]
  parent <- match(node.left/2, node)

  left.child <- match(2 * node, node)
  right.child <- match(node * 2 + 1, node)
  rows <- labels(tree, pretty = pretty)

  xytmp <- FUN2(x = xy$x, y = xy$y, node = node)
  leftptx <- (xytmp$x[2, ] + xytmp$x[1, ])/2
  leftpty <- (xytmp$y[2, ] + xytmp$y[1, ])/2
  rightptx <- (xytmp$x[3, ] + xytmp$x[4, ])/2
  rightpty <- (xytmp$y[3, ] + xytmp$y[4, ])/2
  text(leftptx, leftpty + 0.52 * cxy[2],
       rows[left.child[!is.na(left.child)]], ...)
  text(rightptx, rightpty - 0.52 * cxy[2],
       rows[right.child[!is.na(right.child)]], ...)


  leaves <- if (all) rep(TRUE, nrow(frame))
  else frame$var == "<leaf>"

  if (is.null(frame$yval2))
    stat <- tree$functions$text(yval = frame$yval[leaves], dev = frame$dev[leaves],
                                wt = frame$wt[leaves], ylevel = ylevels, digits = digits,
                                n = frame$n[leaves], use.n = use.n)
  else stat <- tree$functions$text(yval = frame$yval2[leaves, ],
                                   dev = frame$dev[leaves],
                                   wt = frame$wt[leaves], ylevel = ylevels,
                                   digits = digits, n = frame$n[leaves], use.n = use.n)

  oval <- function(middlex, middley, a, b) {
    theta <- seq(0, 2 * pi, pi/30)
    newx <- middlex + a * cos(theta)
    newy <- middley + b * sin(theta)
    polygon(newx, newy, border = TRUE, col = "white")
  }
  rectangle <- function(middlex, middley, a, b) {
    newx <- middlex + c(a, a, -a, -a)
    newy <- middley + c(b, -b, -b, b)
    polygon(newx, newy, border = TRUE, col = "white")
  }

  stat <- gsub(" *$","",stat)
  for (i in parent) oval(xy$x[i], xy$y[i],
                         a = sqrt(2) * 1.2*strwidth(stat[i])/2,
                         b = sqrt(2) * cxy[2])
  child <- match(node[frame$var == "<leaf>"], node)
  for (i in child) rectangle(xy$x[i], xy$y[i],
                             a = 1.2*strwidth(stat[i])/2,
                             b = cxy[2])

  text(xy$x[leaves], xy$y[leaves], stat, adj = 0.5, ...)

  if (!missing(title.))
    title(title., cex = 0.8)
  invisible()
}

