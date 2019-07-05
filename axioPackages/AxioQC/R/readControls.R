readControls <- function( Files = NULL )
{
  RNADat10 <- data.frame()
  RNADat100 <- data.frame()
  DNADat10 <- data.frame()
  DNADat100 <- data.frame()
  for ( i in Files )
  {
    tmpWB <- loadWorkbook( i )
    avSheets <- getSheets( tmpWB )
    sheet10 <- avSheets[grep( "10" , avSheets , fixed = TRUE )]
    sheet100 <- avSheets[grep( "100" , avSheets , fixed = TRUE )]
    sheet10 <- setdiff( sheet10 , sheet100 )
    rna10 <- sheet10[grep( "RNA" , sheet10 , fixed = TRUE )]
    dna10 <- sheet10[grep( "DNA" , sheet10 , fixed = TRUE )]
    rna100 <- sheet100[grep( "RNA" , sheet100 , fixed = TRUE )]
    dna100 <- sheet100[grep( "DNA" , sheet100 , fixed = TRUE )]
    DataFrame <- data.frame( df = c( "RNADat10" , "DNADat10" , "RNADat100" , "DNADat100" ) , Sheet = c( ifelse( length( rna10 ) < 1 , NA , rna10 ) , ifelse( length( dna10 ) < 1 , NA , dna10 ) , ifelse( length( rna100 ) < 1 , NA , rna100 ) , ifelse( length( dna100 ) < 1 , NA , dna100 ) ) )
    dfIdx <- which( !is.na( DataFrame[,2] ) )
    DataFrame <- DataFrame[dfIdx,]
    for ( j in seq( 1 , nrow( DataFrame ) ) )
    {
      eval( parse( text = paste( DataFrame[j,1] , " <- readControlSheet( tmpWB , \"" , DataFrame[j,2] , "\" , controlDat =  " , DataFrame[j,1] , "  )" , sep = "" ) ) )
    }
  }
  if ( length( dna10 ) > 0 )
  {
    DNADat10$QC_Date <- as.Date( as.character( DNADat10$QC_Date ) , format = "%Y-%m-%d" )
    DNADat10$"10ng" = as.numeric( DNADat10$Concentration ) * 1000
  }
  if ( length( dna100 ) > 0 )
  {
    DNADat100$QC_Date <- as.Date( as.character( DNADat100$QC_Date ) , format = "%Y-%m-%d" )
    DNADat100$"100ng" = as.numeric( DNADat100$Concentration ) * 1000
  }
  if ( length( rna10 ) > 0 )
  {
    RNADat10$QC_Date <- as.Date( as.character( RNADat10$QC_Date ) , format = "%Y-%m-%d" )
    RNADat10$"10ng" = as.numeric( RNADat10$Concentration ) * 1000
  }
  if ( length( rna100 ) > 0 )
  {
    RNADat100$QC_Date <- as.Date( as.character( RNADat100$QC_Date ) , format = "%Y-%m-%d" )
    RNADat100$"100ng" = as.numeric( RNADat100$Concentration ) * 1000
  }
  return( list( DNA10 = DNADat10 , DNA100 = DNADat100 , RNA10 = RNADat10 , RNA100 = RNADat100 ) )
}

readControlSheet <- function( controlWB , controlSheet , controlDat = NULL )
{
  tmp <- readWorksheet( controlWB , sheet = controlSheet )
  tmpNames <- names( tmp )
  controlCols <- grep( "barcode" , tmpNames , ignore.case = TRUE )
  controlCols <- c( controlCols , grep( "date" , tmpNames , ignore.case = TRUE ) )
  controlCols <- c( controlCols , grep( "vol" , tmpNames , ignore.case = TRUE ) )
  controlCols <- c( controlCols , grep( "conc" , tmpNames , ignore.case = TRUE ) )
  tmp <- tmp[,controlCols]
  names( tmp ) <- c( "Analysis Barcode" , "QC_Date" , "Volume" , "Concentration" )
  if ( nrow( controlDat ) < 1 )
  {
    controlDat <- tmp
  }
  else
  {
    controlDat <- rbind( controlDat , tmp )
  }
  return( controlDat )
}
