getPostSigTau2 <- function( tau2 , nu , sc , s0 , muSc , sigSc )
{
  logPrior <- dgamma( sc , muSc , sigSc , log = TRUE )
  logLike <- sum( get.dichi2.log( tau2 , nu , sc + s0 ) )
  return( logPrior + logLike )
}
