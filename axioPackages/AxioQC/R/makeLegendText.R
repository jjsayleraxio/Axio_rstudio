makeLegendText <- function( x )
{
  legendText <- vector( length = length( x ) )
  for ( i in seq( 1 , length( x ) ) )
  {
    legendText[i] <- paste( paste( paste( paste( names( x )[i] , ", lower" , sep = "" ) , round( x[[i]][1] , 3 ) ) , "upper" ) , round( x[[i]][2] , 3 ) )
  }
  return( legendText )
}

makeLegendFrame <- function( x )
{
  legendFrame <- data.frame()
  for ( i in seq( 1 , length( x ) ) )
  {
    if ( nrow( legendFrame ) < 1 )
    {
      legendFrame <- data.frame( Date = names( x )[i] , Lower = round( x[[i]][1] , 3 ) , Upper = round( x[[i]][2] , 3 ) )
    }
    else
    {
      legendFrame <- rbind( legendFrame , data.frame( Date = names( x )[i] , Lower = round( x[[i]][1] , 3 ) , Upper = round( x[[i]][2] , 3 ) ) )
    }
  }
  return( legendFrame )
}
