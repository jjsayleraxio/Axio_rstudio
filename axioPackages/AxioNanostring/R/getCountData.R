getCountData <- function( x , multiplexed = FALSE )
{
  if ( multiplexed )
  {
    xNames <- unlist( strsplit( names( x$sampleAttributes )[2] , " " , fixed = TRUE ) )
  }
  else
  {
    xNames <- names( x$sampleAttributes )[2]
  }
  x <- x$codeSummary
  cMat <- matrix( nrow = length( unique( x$Name ) ) , ncol = length( xNames ) )
  rownames( cMat ) <- unique( x$Name )
  colnames( cMat ) <- xNames
  controlIndex <- which( x[,"CodeClass"] %in% c( "Positive" , "Negative" ) )
  controls <- x[controlIndex,"Count"]
  names( controls ) <- x[controlIndex,"Name"]
  x <- x[-controlIndex,]
  x$Group <- as.integer( substr( x$CodeClass , nchar( as.character( x$CodeClass ) ) , nchar( as.character( x$CodeClass ) ) ) )
  for ( i in seq( 1 , length( xNames ) ) )
  {
    indx <- which( x$Group %in% i )
    cMat[c( names( controls ) , as.character( x[indx,"Name"] ) ),i] <- c( controls , x[indx,"Count"] )
  }
  return( cMat )
}
