getPowerElstonPGxNormal <- function( gmeans , sigma2 , p = 0.2 , n = 100 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , type = c( "main" , "interaction" ) )
{
  w <- makeContrastMatrix( gmeans , gaction , type )
  if ( !inherits( gmeans , "matrix" ) )
  {
    gmeans <- as.matrix( gmeans )
  }
  if ( nrow( gmeans ) < 3 )
  {
    gmeans <- t( gmeans )
  }
  gfreq <- c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 )
  n <- matrix( n , nrow = nrow( gmeans ) , ncol = ncol( gmeans ) )
  if ( type == "interaction" )
  {
    if ( length( r < 2 ) )
    {
      r <- c( r , 1 - r )
    }
    for ( i in seq( 1 , 3 ) )
    {
      for ( j in seq( 1 , ncol( gmeans ) ) )
      {
        n[i,j] <- r[j] * gfreq[i] * n[i,j]
      }
    }
  }
  else
  {
    for ( i in seq( 1 , 3 ) )
    {
      n[i,1] <- gfreq[i] * n[i,1]
    }
  }
  Ea <- sum( as.vector( w * gmeans ) ) / sqrt( sigma2 * sum( as.vector( w )^2 / as.vector( n ) ) )
  const <- abs( qnorm( 1 - alpha / 2 ) )
  return( round( pnorm( Ea - const ) + pnorm( -Ea - const ) , 4 ) )
}

getPowerElstonPGxNormalVec <- function( p = 0.2 , gmeans , sigma2 , n = 100 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , type = c( "main" , "interaction" ) )
{
  return( as.vector( unlist( lapply( gmeans , getPowerElstonPGxNormal , sigma2 , p , n , r , gaction , alpha , type ) ) ) )
}

computePowerElstonPGx <- function ( effects , sigma = 1 , p = 0.2 , n = 100 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  stopifnot( inherits( effects , "list" ) )
  powerMatrix <- matrix( unlist( sapply( p , getPowerElstonPGxNormalVec , effects , sigma^2 , n , r , gaction , alpha , type ) ) , byrow = TRUE , nrow = length( p ) , ncol = length( effects ) )
  sampleMatrix <- matrix( unlist( sapply( p , getSampleElstonPGxNormalVec , effects , sigma^2 , r , gaction , alpha , beta , type ) ) , byrow = TRUE , nrow = length( p ) , ncol = length( effects ) )
  return( new( "powerObject" , effects = effects,
          MAF = p,
          power = powerMatrix,
          sampleSize = sampleMatrix,
          SD = sigma,
          N = n,
          alpha = alpha,
          beta = beta,
          proportion = r,
          type = type ,
          stype = "Normal" ) )
}
