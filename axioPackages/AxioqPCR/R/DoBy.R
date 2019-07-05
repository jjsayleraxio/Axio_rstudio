DoBy <- function( response = NULL , group = NULL , data = NULL , FUN = "mean" , ... )
{
  rVec <- vector()
  z <- vector()
  rMatNames <- vector()
  classVec <- vector()
  for ( i in seq( length( group ) ) )
  {
    classVec[i] <- class( data[,group[i]] )
  }
  for ( i in seq( length( FUN ) ) )
  {
    for ( j in seq( length( response ) ) )
    {
      z <- tapply( data[,response[j]] , data[,group] , eval( FUN[i] ) , ... )
      rVec <- cbind( rVec , as.vector( z ) )
      rMatNames <- rbind( rMatNames , paste( response[j] , FUN[i] , sep = "." ) )
    }
  }
  d <- dimnames( z )
  rMat <- cbind( expand.grid( d ) , rVec )
  if ( length( group ) > 1 )
  {
    names( rMat ) <- c( names( d ) , rMatNames )
  }
  else
  {
    names( rMat ) <- c( group , rMatNames )
  }
  for ( i in seq( length( group ) ) )
  {
    if ( classVec[i] == "POSIXct" )
    {
      rMat[,group[i]] <- as.POSIXct( rMat[,group[i]] )
    }
    else if ( classVec[i] == "Date" )
    {
      rMat[,group[i]] <- as.Date( rMat[,group[i]] )
    }
    else
    {
      rMat[,group[i]] <- as( as.character( rMat[,group[i]] ) , classVec[i] )
    }
  }
  return( rMat )
}
