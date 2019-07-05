readFails <- function( Files = NULL , failsNames = c( "Date" , "Batch ID" , "Job ID" , "Biomek A#" , "control barcode" , "high/low control?" , "concentration" , "vol" , "well" ) )
{
  rnaFails <- data.frame()
  dnaFails <- data.frame()
  for ( i in Files )
  {
    tmpWB <- loadWorkbook( i )
    failtmp <- readWorksheet( tmpWB , sheet = "Table" )
    failtmp <- failDataFrame( failtmp , failsNames )
    if ( nrow( rnaFails ) < 1 )
    {
      rnaFails <- failtmp[grep( "US" , failtmp$"control barcode" , fixed = TRUE ),]
    }
    else
    {
      rnaFails <- rbind( rnaFails , failtmp[grep( "US" , failtmp$"control barcode" , fixed = TRUE ),] )
    }
    if ( nrow( dnaFails ) < 1 )
    {
      dnaFails <- failtmp[grep( "DS" , failtmp$"control barcode" , fixed = TRUE ),]
    }
    else
    {
      dnaFails <- rbind( dnaFails , failtmp[grep( "DS" , failtmp$"control barcode" , fixed = TRUE ),] )
    }
  }
  return( list( DNA10 = data.frame( "Analysis Barcode" = dnaFails$"control barcode" , "QC_Date" = dnaFails$Date , Volume = dnaFails$vol , Concentration = dnaFails$concentration / 1000 , "10ng" = dnaFails$concentration , check.names = FALSE )[which( dnaFails$concentration < 30 ),] , DNA100 = data.frame( "Analysis Barcode" = dnaFails$"control barcode" , "QC_Date" = dnaFails$Date , Volume = dnaFails$vol , Concentration = dnaFails$concentration / 1000 , "100ng" = dnaFails$concentration , check.names = FALSE )[which( dnaFails$concentration >= 30 ),] , RNA10 = data.frame( "Analysis Barcode" = rnaFails$"control barcode" , "QC_Date" = rnaFails$Date , Volume = rnaFails$vol , Concentration = rnaFails$concentration / 1000 , "10ng" = rnaFails$concentration , check.names = FALSE )[which( rnaFails$concentration < 30 ),] , RNA100 = data.frame( "Analysis Barcode" = rnaFails$"control barcode" , "QC_Date" = rnaFails$Date , Volume = rnaFails$vol , Concentration = rnaFails$concentration / 1000 , "100ng" = rnaFails$concentration , check.names = FALSE )[which( rnaFails$concentration >= 30 ),] ) )
}
