getAnnotationData <- function( x , multiplexed = FALSE )
{
  if ( multiplexed )
  {
    xNames <- unlist( strsplit( names( x$sampleAttributes )[2] , " " , fixed = TRUE ) )
  }
  else
  {
    xNames <- names( x$sampleAttributes )[2]
  }
  a <- x$sampleAttributes
  b <- x$laneAttributes
  names( a ) <- names( b ) <- c( "A" , "B" )
  x <- rbind( a , b )
  cMat <- data.frame( matrix( nrow = length( xNames ) , ncol = nrow( x ) ) )
  rownames( cMat ) <- xNames
  colnames( cMat ) <- x[,1]
  for ( i in seq( 1 , length( xNames ) ) )
  {
    cMat[i,] <- x[,2]
  }
  cMat$Date <- as.Date( cMat$Date , format = "%Y%m%d" )
  return( cMat )
}
