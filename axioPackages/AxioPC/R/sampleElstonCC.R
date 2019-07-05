
getSampleElstonPGxCC <- function( gprobs , p = 0.2 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  w <- makeContrastMatrix( gprobs , gaction , type )
  if ( !inherits( gprobs , "matrix" ) )
  {
    gprobs <- as.matrix( gprobs )
  }
  gprobs <- t( gprobs )
  if ( nrow( gprobs ) < 3 )
  {
    gprobs <- t( gprobs )
  }
  n <- matrix( NA , nrow = nrow( gprobs ) , ncol = ncol( gprobs ) )
  if ( type == "interaction" )
  {
    for ( i in seq( 1 , ncol( n ) ) )
    {
      n[,i] <- r[i] * c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 )
    }
  }
  else
  {
    n <- matrix( c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 ) , nrow = nrow( gprobs ) )
  }
  return( round( ( ( qnorm( 1 - alpha / 2 ) + qnorm( beta ) ) / sum( as.vector( w * gprobs ) ) )^2 * ( sum( as.vector( w )^2 * as.vector( gprobs ) * ( 1 - as.vector( gprobs ) ) / as.vector( n ) ) ) ) )
}

getSampleElstonPGxCCVec <- function( p = 0.2 , gprobs , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  return( as.vector( lapply( gprobs , getSampleElstonPGxCC , p , r , gaction , alpha , beta , type ) ) )
}
