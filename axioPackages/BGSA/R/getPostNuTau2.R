getPostNuTau2 <-
function(tau2, nu, sc, muNu, sigNu) {
	logPrior = dgamma(nu, muNu, sigNu, log=TRUE)
	logLike = sum(get.dichi2.log(tau2, nu, sc))
	return(logPrior+logLike)
}
