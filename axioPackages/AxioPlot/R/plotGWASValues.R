plotGWASValues <- function( position , chromosome , values , geneIDs , sort = TRUE , threshold = NULL , cutoff = NULL , suggestive = NULL , threshLabels = FALSE , geneLabels = FALSE , yLabel = NULL , mainTitle = NULL , ... )
{
  if ( !is.null( cutoff ) )
  {
    coff <- values > -log10( cutoff ) & !is.na( values )
    position <- position[coff]
    chromosome <- chromosome[coff]
    values <- values[coff]
    geneIDs <- geneIDs[coff]
  }
  if ( sort )
  {
    ordered <- order(chromosome,position)
    position <- position[ordered]
    chromosome <- chromosome[ordered]
    values <- values[ordered]
    geneIDs <- geneIDs[ordered]
  }
  NewPos <- cumsum( as.numeric( position ) )
  
  xTicks <- unique( chromosome , na.rm = TRUE )
  xTicks <- xTicks[order(xTicks)]
  xNames <- xTicks
  names( xTicks ) <- xTicks
  names( xNames ) <- xTicks
  color <- 0
  Color <- vector( length = length( ordered ) )
  
  for ( i in seq_along( xTicks ) )
  {
    Color[chromosome == xTicks[i]] <- color
    xTicks[i] <- max( NewPos[chromosome == xTicks[i]] )
    xNames[i] <- ( min( NewPos[chromosome == xNames[i]] , na.rm = TRUE ) + max( NewPos[chromosome == xNames[i]] ) ) / 2
    color <- abs( color - 1 )
  }
  xTicks <- c( min( NewPos ) , xTicks )
  Color <- Color + 1
  
  if ( is.null( mainTitle ) )
  {
    mainTitle <- "GWAS Plot"
  }
  xLabel <- "Chromosome Map Position (MB)"
  if ( is.null( yLabel ) )
  {
    yLabel <- expression(-Log[10]~~group("(",P-Value,")"))
  }
  if ( is.null( threshold ) )
  {
    yLim <- c( min(values,na.rm=TRUE)-0.1 , max(values)+0.1 )
  }
  else
  {
    yLim <- c( min(values,na.rm=TRUE)-0.1 , c( max(values)+0.1 , -log10( threshold )+0.1 )[which.max( c( max(values)+0.1 , -log10( threshold )+0.1 ) )] )
  }
  plot( NewPos , values , pch=20 , col = Color , main = mainTitle , xlab = xLabel , ylab = yLabel , ylim = yLim , axes = FALSE , type = "l" , ... )
  axis( 2 )
  axis( 1 , xTicks , labels = FALSE )
  axis( 1 , xNames , names( xNames ) , tick = FALSE )
  if ( !is.null( threshold ) )
  {
    sigGenes <- values > -log10( threshold )
    if ( geneLabels )
    {
      if ( any( sigGenes ) )
      {
        text( NewPos[sigGenes] , pValue[sigGenes] , geneIDs[sigGenes] , pos = 4 , offset = .3 )
      }
      else
      {
        warning( "No Genes Significant at p < " , threshold )
      }
    }
    abline( h = -log10( threshold ) , col = "red" , lty = 3 )
    if ( threshLabels )
    {
      text( mean( xTicks[2:3] ) , -log10( threshold ) + 0.15 , paste( "(Genome-wide significant, p < " , threshold , ")" , sep = "" ) , pos = 2 , col = "blue" )
    }
  }
  if ( !is.null( suggestive ) )
  {
    abline( h = -log10( suggestive ) , col = "red" , lty = 2 )
    if ( threshLabels )
    {
      text( mean( xTicks[2:3] ) , -log10( suggestive ) + 0.15 , paste( "(Genome-wide suggestive, p < " , suggestive , ")" , sep = "" ) , pos = 2 , col = "blue" )
    }
  }
}
