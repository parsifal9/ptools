match.prior <- function(prior = c("NG", "R", "DEG", "NEG", "SCAD")) {
	prior <- match.arg(prior)

	cexpfn <- switch(prior,
		NG   = HGcexp.ng,
		R    = HGcexp.r,
		DEG  = HGcexp.deg,
		NEG  = HGcexp.neg,
		SCAD = HGcexp.scad)

	cexpfn
}
