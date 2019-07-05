plotEQTL <- function( eqtlLocCounts = NULL , filename = "" , lodCutOff = NULL , sort = TRUE , type = c( "points" , "hexbin" ) , threshold = NULL , cutoff = NULL , suggestive = NULL , threshLabels = FALSE , geneLabels = FALSE , mainTitle = NULL , ... )
{
	if ( is.null( eqtlLocCounts ) )
  {
    if ( identical( filename , "" ) )
    {
      stop( "File name must not be empty when eqtlLocCounts is NULL" )
    }
    gawkCMD <- paste( paste( "gawk \'$5>23.03 {print $1, $2, $3, $4, $5}\'" , filename ) , "> .tmp.txt" )
	  system( gawkCMD )
	  eqtl <- read.table( ".tmp.txt" , header = TRUE , as.is = TRUE )
	  system( "rm .tmp.txt" )
	  names( eqtl ) <- c( "gene" , "chr" , "marker" , "pos" , "LLR" )
	  eqtl$LOD <- eqtl$LLR / 4.6
	  if ( !is.null( lodCutOff ) )
	  {
	    eqtl <- eqtl[eqtl$LOD > lodCutOff,]
	  }
	  eqtlLocCounts <- as.data.frame( table( eqtl$marker ) )

	  eqtlLocCounts <- merge.data.frame( eqtl[!duplicated( eqtl$marker ) , c( "chr" , "pos" , "marker" )] , eqtlLocCounts , by.x = 3 , by.y = 1 )
  }

	eqtlLocCounts <- eqtlLocCounts[order( eqtlLocCounts$chr , eqtlLocCounts$pos ),]

	plotGWASValues( position = eqtlLocCounts$pos , chromosome = eqtlLocCounts$chr , values = eqtlLocCounts$Freq , sort = FALSE , yLabel = ifelse( is.null( lodCutOff ) , "eQTL Counts" , paste( "eQTL Counts >" , lodCutOff ) ) , mainTitle = "eQTL by Location" , ... )
}
