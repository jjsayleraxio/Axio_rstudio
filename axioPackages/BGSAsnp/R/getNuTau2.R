getNuTau2 <- function( tau2 , nu , sc , muNu , sigNu )
{
  w <- 1
  m <- 100
  z <- getPostNuTau2( tau2 , nu , sc , muNu , sigNu ) - rexp( 1 )
  # Stepping out to obtain the [L, R] range
  u <- runif (1 )
  L <- nu - w * u
  R <- L + w
  v <- runif( 1 )
  J <- floor( m * v )
  K <- ( m - 1 ) - J
  L <- max( 0 , L )
  while (J > 0 && L > 0 && z < getPostNuTau2( tau2 , L , sc , muNu , sigNu ) )
  {
    L <- L - w  
    L <- max( 0 , L )    
    J <- J - 1
  }
  while ( K > 0 && z < getPostNuTau2( tau2 , nu , R , muNu , sigNu ) )
  {
    R <- R + w
    K <- K - 1
  }
  # Shrinkage to obtain a sample
  u <- runif( 1 )
  newParam <- L + u * (R - L )
  while ( z > getPostNuTau2( tau2 , newParam , sc , muNu , sigNu ) )
  {
    if ( newParam < nu )
    {
      L <- newParam
    }
    else
    {
      R <- newParam
    }
    u <- runif( 1 )
    newParam <- L + u * ( R - L )
  }
  return( newParam )
}
