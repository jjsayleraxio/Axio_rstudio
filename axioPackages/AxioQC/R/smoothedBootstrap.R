smoothedBootstrap <- function( x , n )
{
  rMat <- matrix( NA , nrow = length( x ) , ncol = n )
  scale <- ifelse( length( x ) < 2 , 1 , sd( x ) * sqrt( 2.0 ) / ( length( x ) )^0.75 )
  for ( i in seq( 1 , n ) )
  {
    rMat[,i] <- sample( x , length( x ) , replace = TRUE ) + scale * rnorm( length( x ) )
  }
  return( rMat )
}
