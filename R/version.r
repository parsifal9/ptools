HGmultc.version <- function() {
	wh <- .path.package("HGmultc")
	tmp <- readLines(file.path(wh, "DESCRIPTION"))
	grep("^Version:", tmp, value = TRUE)
}
