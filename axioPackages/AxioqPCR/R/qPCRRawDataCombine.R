qPCRRawDataCombine <- function( sdsDirectory = getwd() , outputFileName = "" , qPCRPlatform = c( "STD" , "TLDA" ) , fileList = NULL , excel = FALSE , thresholdLimit = 0.2 , undetermined = NULL )
{
  if ( identical( R.home() , sdsDirectory ) )
  {
    error.msg <- paste( "You seem to have set your working directory to your R program directory.",
        "To change working directory to the folder containing the raw qPCR data saved in .txt format,\n" ,
        "specify a different directory for the sdsDirectory argument." , sep = "" )
    stop( error.msg )
  }
#  if ( qPCRPlatform != "STD" & qPCRPlatform != "TLDA" | qPCRPlatform == "" )
#  {
#    error.msg <- "qPCR platform must be either \"STD\" or \"TLDA\".\n"
#    stop( error.msg )
#  }
  if ( outputFileName == "" )
  {
    error.msg <- "You need to enter a file name for the final file."
    stop( error.msg )
  }
  qPCRPlatform <- toupper( match.arg( qPCRPlatform ) )
  
  ## Find all the .txt files
  if ( is.null( fileList ) )
  {
    dataList <- file.path( sdsDirectory , list.files( sdsDirectory , pattern = ".txt" ) )
  }
  else
  {
    dataList <- file.path( sdsDirectory , fileList )
  }
  
  ## combine all raw qPCR data into a .csv file
  output <- data.frame()
  for ( i in seq( 1 , length( dataList ) ) )
  {
    cat( "Starting file:" , dataList[i] , "..." )
    temp <- readSDStxt( dataList[i] , qPCRPlatform , i , undetermined )
    output <- rbind( output , temp )
    cat( " Done!\n" )
  }
  
  baselineType <- unique( output$"Baseline Type" )
  thresholdType <- unique( output$"Threshold Type" )
  threshold <- unique( output$Threshold )
  # Stop if qPCR analysis settings were wrong
  if ( baselineType != "Automatic" & ( thresholdType != "Manual" | threshold != thresholdLimit ) )
  {
    error.msg <- "Analysis settings for SDS file were not set correct.  Re-set analysis settings and re-export data in .txt format."
    stop( error.msg )
    break
  }
  
  output <- output[,-which( names( output ) %in% c( "Baseline Type" , "Baseline Start" , "Baseline Stop" , "Threshold Type" , "Threshold" ) )]
  
  if ( excel )
  {
    sdsSetting <- matrix( c( baselineType ,thresholdType ,threshold ) , ncol = 1 )
    rownames( sdsSetting ) <- c( "Baseline Type" , "Threshold Type" , "Threshold" )
    require( XLConnect )
    outputWorkbook <- loadWorkbook( file.path( pathRoot , outputFileName ) , create = TRUE )
    createSheet( outputWorkbook , name = "CTs" )
    writeWorksheet( object = outputWorkbook , data = output , sheet = "CTs" , header = TRUE )
    createSheet( outputWorkbook , name = "SDS Setting" )
    writeWorksheet( object = outputWorkbook , data = cbind( rownames( sdsSetting) , sdsSetting ) , sheet = "SDS Setting" , header = FALSE )
    saveWorkbook( outputWorkbook )
  }
  else
  {
    sdsSetting <- paste( "Baseline Type: " , baselineType , "\nThreshold Type: " , thresholdType , "\nThreshold: " , threshold , sep = "" )
    write.table( sdsSetting , file = file.path( pathRoot , "SDS_Setting.csv" ) , row.names = FALSE , quote = FALSE , col.names = FALSE )
    write.csv( output , file = file.path( pathRoot , outputFileName ) , row.names = FALSE )
  }
  invisible()
}
