failDataFrame <- function( x , failsNames )
{
  x <- x[-which( apply( apply( x , 1 , is.na ) , 2 , all ) ),]
  failsList <- vector( mode = "list" )
  for ( i in failsNames )
  {
    failsList[[i]] <- x[which( x[,1] %in% i ),-1]
  }
  x <- data.frame( lapply( failsList , function( x ) matrix( unlist( x ) , ncol = 1 ) ) )
  x <- x[-which( apply( apply( x , 1 , is.na ) , 2 , all ) ),]
  names( x ) <- failsNames
  x <- x[-which( is.na( x$Date ) ),]
  tmp <- as.Date( x$Date , format = "%Y-%m-%d %H:%M:%S" )
  tmp[which( is.na( tmp ) )] <- as.Date( as.character( x$Date[which( is.na( tmp ) )] ) , format = "%m/%d/%Y" )
  x$Date <- tmp
  x$concentration <- as.numeric( as.character( x$concentration ) )
  x$vol <- as.numeric( as.character( x$vol ) )
  x$"control barcode" <- toupper( as.character( x$"control barcode" ) )
  x$"well" <- toupper( as.character( x$"well" ) )
  x$"high/low control?" <- factor( tolower( as.character( x$"high/low control?" ) ) )
  x <- x[-which( is.na( x$Date ) | is.na( x$concentration ) | x$concentration == 0.0 | is.na( x$vol ) ),]
  return( x )
}
