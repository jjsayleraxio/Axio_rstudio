
getSampleElstonPGxNormal <- function( gmeans , sigma2 , p = 0.2 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  w <- makeContrastMatrix( gmeans , gaction , type )
  if ( !inherits( gmeans , "matrix" ) )
  {
    gmeans <- as.matrix( gmeans )
  }
  gmeans <- t( gmeans )
  if ( nrow( gmeans ) < 3 )
  {
    gmeans <- t( gmeans )
  }
  n <- matrix( NA , nrow = nrow( gmeans ) , ncol = ncol( gmeans ) )
  if ( type == "interaction" )
  {
    if ( length( r ) < 2 )
    {
      r <- c( r , 1 - r )
    }
    for ( i in seq( 1 , ncol( n ) ) )
    {
      n[,i] <- r[i] * c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 )
    }
  }
  else
  {
    n <- matrix( c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 ) , nrow = nrow( gmeans ) )
  }
  return( round( ( ( qnorm( 1 - alpha / 2 ) + qnorm( beta ) ) / sum( as.vector( w * gmeans ) ) )^2 * ( sigma2 * sum( as.vector( w )^2 / as.vector( n ) ) ) ) )
}

getSampleElstonPGxNormalVec <- function( p = 0.2 , gmeans , sigma2 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  return( as.vector( lapply( gmeans , getSampleElstonPGxNormal , sigma2 , p , r , gaction , alpha , beta , type ) ) )
}
