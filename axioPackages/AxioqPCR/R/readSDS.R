readSDS <- function( dataFile = "" , qPCRPlatform = "STD" , plate = NULL )
{
  sdsColClasses <- c( rep( "character" , 3 ) , rep( "NULL" , 2 ) , rep( "character" , 12 ) , rep( "NULL" , 13 ) )
  sdsVals <- read.table( dataFile , nrows = 4 , stringsAsFactor = FALSE , sep = "\t" , header = TRUE , fill = TRUE , colClasses = c( "NULL" , "character" , "NULL" ) )[c( 1 , 2 , 4 ),]
  names( sdsVals ) <- c( "FileName" , "PlateID" , "RunDateTime" )
  temp <- read.table( dataFile , skip = 10 , header = TRUE , fill = TRUE , sep = "\t" , stringsAsFactor = FALSE , colClasses = sdsColClasses , na.strings = "Undetermined" )
  ## Stop if non-SDS text file is in the folder.
  if ( colnames( temp )[1] != "Well" )
  {
    errorMsg <- paste( dataFile , "does not contain the right ABI SDS format.\n  Remove the file from the folder." )
    stop( errorMsg )
    break
  }
  rowsKeepMk <- which( !( temp$Well %in% c( "" , "Slope" , "Y-Intercept" , "R^2" , "NAP" , "Well" ) | temp$Detector.Name %in% "" ) )
  temp <- temp[rowsKeepMk,]
  sdsColClasses <- c( "integer" , rep( "character" , 2 ) , rep( "numeric" , 7 ) , rep( "character" , 4 ) , "numeric" )
  for ( j in seq( 1 , ncol( temp ) ) )
  {
    tfunc <- paste( "temp[,j] <- as." ,  sdsColClasses[j] , "( temp[,j] )" , sep = "" )
    eval( parse( text = tfunc ) )
  }
  if ( qPCRPlatform == "STD" )
  {
    temp$File.Name <- sdsVals["FileName"]
    temp$Plate.ID <- sdsVals["PlateID"]
    temp$Run.DateTime <- sdsVals["RunDateTime"]
  }
  else
  {
    temp$File.Name <- sdsVals["FileName"]
    temp$Plate.ID <- sdsVals["PlateID"]
    temp$Run.DateTime <- sdsVals["RunDateTime"]
    temp$Plate.No <- plate
  }
  return( temp )
}

readSDStxt <- function( dataFile = "" , qPCRPlatform = "STD" , plate = NULL , undetermined = NULL )
{
  oobNames <- c( NA , "" , "True" , "False" , "SDS 2.3" , "Assay Type" , "Operator" , "ThermalCycleParams" , "Slope" , "Y-Intercept" , "R^2" , "NAP" )
  indx <- which( c( rep( TRUE , 3 ) , rep( FALSE , 2 ) , rep( TRUE , 12 ) , rep( FALSE , 13 ) ) )
  dFile <- file( dataFile , open = "r" )
  sampInfo <- FALSE
  temp <- data.frame()
  tempStarted <- FALSE
  if ( qPCRPlatform == "STD" )
  {
    sdsVals <- vector( length = 3 )
    names( sdsVals ) <- c( "Filename" , "PlateID" , "Run DateTime" )
  }
  else
  {
    sdsVals <- vector( length = 4 )
    names( sdsVals ) <- c( "Filename" , "PlateID" , "Run DateTime" , "Plate" )
  }
  while( length( tmp <- readLines( dFile , n = 1 , warn = FALSE ) ) > 0 )
  {
    tmp <- as.data.frame( matrix( unlist( strsplit( tmp , "\t" , fixed = TRUE ) ) , nrow = 1 ) , stringsAsFactors = FALSE )
    if ( ncol( tmp ) > 0 )
    {
      if ( !( tmp[,1] %in% c( oobNames , "Filename" , "PlateID" , "Run DateTime" , "Plate" ) ) )
      {
        if ( tmp[,1] %in% "Sample Information" )
        {
          sampInfo <- TRUE
        }
        else
        {
          tmp <- tmp[,indx]
          if ( sampInfo )
          {
            dHeader <- as.vector( tmp )
            oobNames <- c( oobNames , "Well" )
            sampInfo <- FALSE
            if ( dHeader[1] != "Well" )
            {
              errorMsg <- paste( dataFile , "does not contain the right ABI SDS format.\n  Remove the file from the folder." )
              stop( errorMsg )
              break
            }
          }
          else
          {
            if ( tempStarted )
            {
              temp <- rbind( temp , tmp )
            }
            else
            {
              temp <- tmp
              tempStarted <- TRUE
            }
          }
        }
      }
      else if( tmp[,1] %in% c( "Filename" , "PlateID" , "Run DateTime" ) )
      {
        sdsVals[tmp[,1]] <- tmp[,2]
      }
      else if( qPCRPlatform == "TLDA" && tmp[,1] %in% "Plate" )
      {
        sdsVals[tmp[,1]] <- tmp[,2]
      }
    }
  }
  close( dFile )
  names( temp ) <- dHeader
  sdsColClasses <- c( "integer" , rep( "character" , 2 ) , rep( "numeric" , 7 ) , rep( "character" , 4 ) , "numeric" )
  for ( j in seq( 1 , ncol( temp ) ) )
  {
    if ( any( temp[,j] %in% "Undetermined" ) )
    {
      temp[which( temp[,j] %in% "Undetermined" ),j] <- ifelse( is.null( undetermined ) , NA , undetermined )
    }
    tfunc <- paste( "temp[,j] <- as." ,  sdsColClasses[j] , "( temp[,j] )" , sep = "" )
    eval( parse( text = tfunc ) )
  }
  if ( qPCRPlatform == "STD" )
  {
    temp$Filename <- sdsVals["Filename"]
    temp$PlateID <- sdsVals["PlateID"]
    temp$"Run DateTime" <- sdsVals["Run DateTime"]
  }
  else
  {
    temp$Filename <- sdsVals["Filename"]
    temp$PlateID <- sdsVals["PlateID"]
    temp$"Run DateTime" <- sdsVals["Run DateTime"]
    temp$"Plate No" <- plate
  }
  return( temp )
}
