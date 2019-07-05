# pg = prob of risk class... e.g. 1-(1-p)^2 for dominant MOI with allele freq p
# pd = prob of disease in population
# OR = odds ratio (may be vector)
# N = sample size (may be vector)
# alpha = significance 

computePowerCohort <- function( pg , pd , OR , N , alpha )
{
  p1 <- ( ( 1 + ( OR - 1 ) * ( pg + pd ) ) - sqrt( ( 1 + ( OR - 1 ) * ( pg + pd ) )^2 - 4 * pd * pg * ( OR - 1 ) * OR ) ) / ( 2 * pg * ( OR - 1 ) )
  p0 <- p1 / ( ( 1 - p1 ) * OR + p1 )
  ES <- sqrt( ( ( p1 - pd ) * pg )^2 / ( pd * pg ) + ( ( p0 - pd ) * ( 1 - pg ) )^2 / ( pd * ( 1 - pg ) ) + ( ( pd - p1 ) * pg )^2 / ( ( 1 - pd ) * pg ) + ( ( pd - p0 ) * ( 1 - pg ) )^2 / ( ( 1 - pd ) * ( 1 - pg ) ) )
  PowOr <- OR
  for ( n in N )
  {
    pow <- pwr.chisq.test( w = ES , N = n , df = 1 , sig.level = alpha )$power
    PowOr <- cbind( PowOr , pow )
  }
  PowOr <- as.data.frame( PowOr )
  names( PowOr) <- c( "OR" , paste( "pow" , N , sep="." ) )
  return( PowOr )
}

computeES <- function( x , pNull )
{
  return( ES.w1( pNull , x ) )
}

# N = (n.case, n.control), may be k-by-2 matrix
computePowerCC <- function( pg , pd , OR , N , alpha )
{
  p1 <- ( ( 1 + ( OR - 1 ) * ( pg + pd ) ) - sqrt( ( 1 + ( OR - 1 ) * ( pg + pd ) )^2 - 4 * pd * pg * ( OR - 1 ) * OR ) ) / ( 2 * pg * ( OR - 1 ) )
  p0 <- p1 / ( ( 1 - p1 ) * OR + p1 )
  pgD <- p1 * pg / pd
  pgND <- ( 1 - p1 ) * pg / ( 1 - pd )
  ES <- 2 * asin( sqrt( pgD ) ) - 2 * asin( sqrt( pgND ) )
  PowOR <- OR
  if ( !is.null( dim( N ) ) )
  {
    for( k in 1:dim( N )[1] )
    {
      nCase <- N[k,1]
      nControl <- N[k,2]

      pow <- pwr.2p2n.test( ES , n1 = nCase , n2 = nControl , sig.level = alpha , alternative = "two.sided" )$power
      PowOR <- cbind( PowOR , pow )
    }
    PowOR <- as.data.frame( PowOR )
    names( PowOR) <- c( "OR ", paste( "pow" , N[,1] , N[,2] , sep="." ) )
  }
  else
  {
    pow <- pwr.2p2n.test( ES , n1 = N[1] , n2 = N[2] , sig.level = alpha , alternative = "two.sided" )$power
    PowOR <- cbind( PowOR , pow )
    PowOR <- as.data.frame( PowOR )
    names( PowOR ) <- c( "OR" , paste( "pow" , N[1] , N[2] , sep="." ) )
  }
  return( PowOR )
}
