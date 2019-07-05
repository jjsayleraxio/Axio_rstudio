sampleKellyPowerNormal <- function( n , gmeans , sigma2 , p , simSamps , r , gaction , alpha , type )
{
  return( getSampleKellyPGxNormal( gmeans , sigma2 , p , n , simSamps , r , gaction , alpha , type ) )
}

getSampleKellyPGxNormal <- function( gmeans , sigma2 , p = 0.2 , r = 0.5 , simSamps = 1000 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  low <- 30
  high <- 200
  ns <- vector()
  while( ( high - low ) > 1 )
  {
    ns <- sample( seq( low , high ) , min( 30 , high - low ) )
    nVec <- as.vector( unlist( lapply( ns , sampleKellyPowerNormal , gmeans , sigma2 , p , simSamps , r , gaction , alpha , type ) ) )
    nVecTmp <- nVec - beta
    print( nVecTmp )
    nVecTmp[which( nVecTmp > 0 )] <- NA
    low <- ns[which.max( nVecTmp )]
    nVecTmp <- nVec - beta
    nVecTmp[which( nVecTmp <= 0 )] <- NA
    high <- ns[which.min( nVecTmp )]
    if ( identical( high , integer( 0 ) ) )
    {
      high <- max( ns ) * 10
    }
    if ( identical( low , integer( 0 ) ) )
    {
      low <- max( trunc( min( ns ) / 10 ) , 1 )
    }
    print( high )
    print( low )
    print( ( high - low ) > 1 )
  }
  return( ns[which.min( ns )] )
}

getSampleKellyPGxNormalVec <- function( p = 0.2 , gmeans , sigma2 , r = 0.5 , simSamps = 1000 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  return( as.vector( lapply( gmeans , getSampleKellyPGxNormal , sigma2 , p , r , simSamps , gaction , alpha , beta , type ) ) )
}
