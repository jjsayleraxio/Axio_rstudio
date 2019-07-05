#######################
# This is the power function for continuous outcome to detect interaction term b3
# Y = b0 + b1 Z + b2 G + b3 G * Z
# two arm: intervention and placebo
# N is total sample size
# r1 is proportion in intervention
# ey0, ey1 is the observed outcome mean in the control arm and in the intervention arm
# sy0, sy1 is the observed outcome SD in the control arm and in the intervention arm
# p is the minor allele frequency
#######################

getPower <- function( sy1 , sy0 , p , b2 , b3 , alpha , N , r1 )
{
  r0 <- 1 - r1
  N1 <- N * r1
  N0 <- N * r0
  var1 <- sy1^2
  var0 <- sy0^2
  sigma2 <- r1 * var1 + r0 * var0 - 2 * p * ( 1 - p ) * ( b2^2 + 2 * r1 * b2 * b3 + r1 * b3^2 )
  coef <- ( 1 / N1 + 1 / N0 ) / ( 2 * p * ( 1 - p ) )
  var <- coef * sigma2
  sigma <- sqrt( var )
  const <- abs( qnorm ( 1 - alpha / 2 ) )
  return( pnorm( -const - b3 / sigma ) + pnorm( -const + b3 / sigma ) )
}

#######################

getSample <- function( sy1 , sy0 , p , b2 , b3 , alpha , power , r1 )
{
  r0 <- 1 - r1
  var1 <- sy1^2
  var0 <- sy0^2
  sigma2 <- r1 * var1 + r0 * var0 - 2 * p * ( 1 - p ) * ( b2^2 + 2 * r1 * b2 * b3 + r1 * b3^2 )
  const <- abs( qnorm( 1 - alpha / 2 ) )
  b <- qnorm( power )
  return( ( b + const ) * ( b + const ) * sigma2 / ( r1 * r0 * b3^2 * 2 * p * ( 1 - p ) ) )
}
