
getPowerKellyPGxCC <- function( gprobs , p = 0.2 , n = 100 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , type = c( "main" , "interaction" ) )
{
  w <- makeContrastMatrix( gprobs , gaction , type )
  if ( !inherits( gprobs , "matrix" ) )
  {
    gmeans <- as.matrix( gprobs )
  }
  if ( nrow( gprobs ) < 3 )
  {
    gprobs <- t( gprobs )
  }
  gfreq <- c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 )
  n <- matrix( n , nrow = nrow( gprobs ) , ncol = ncol( gprobs ) )
  if ( type == "interaction" )
  {
    if ( length( r < 2 ) )
    {
      r <- c( r , 1 - r )
    }
    for ( i in seq( 1 , 3 ) )
    {
      for ( j in seq( 1 , ncol( gprobs ) ) )
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
  Ea <- sum( as.vector( w * gprobs ) ) / sqrt( sum( as.vector( w )^2 * as.vector( gprobs ) * ( 1 - as.vector( gprobs ) ) / as.vector( n ) ) )
  const <- abs( qnorm( 1 - alpha / 2 ) )
  return( round( pnorm( Ea - const ) + pnorm( -Ea - const ) , 4 ) )
}

getPowerKellyPGxCCVec <- function( p = 0.2 , gprobs , n = 100 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , type = c( "main" , "interaction" ) )
{
  return( as.vector( lapply( gprobs , getPowerKellyPGxCC , p , n , r , gaction , alpha , type ) ) )
}

computePowerKellyCC <- function ( effects , p = 0.2 , n = 100 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  stopifnot( inherits( effects , "list" ) )
  powerMatrix <- matrix( nrow = length( p ) , ncol = length( effects ) )
  sampleMatrix <- matrix( nrow = length( p ) , ncol = length( effects ) )
  powerMatrix <- sapply( p , getPowerKellyPGxCCVec , effects , n , r , gaction , alpha , type )
  sampleMatrix <- sapply( p , getSampleKellyPGxCCVec , effects , r , gaction , alpha , beta , type )
  return( new( "powerObject" , effects = effects ,
          MAF = p ,
          power = powerMatrix ,
          sampleSize = sampleMatrix ,
          SD = 0 ,
          N = n ,
          alpha = alpha ,
          beta = beta ,
          proportion = r ,
          type = type ,
          stype = "Case Control" ) )
}
