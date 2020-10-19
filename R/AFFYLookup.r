AFFYLookup <- function(probe.list, unigene = TRUE) {
	if (!exists("AffyU95genes"))
	  	data(AffyU95genes)

	probes <- "Probe.set.Name"
	if (unigene)
		annot <- "Unigene.Title"
	else
		annot <- "Affy.Descriptions"

	n <- length(AffyU95genes[, annot])
	tt <- match(make.names(probe.list),
		make.names(AffyU95genes[, probes]),
		nomatch = n + 1)

	xx <- c(AffyU95genes[, annot], "Not Found")[tt]
	names(xx) <- make.names(probe.list)
	return(xx)
}
