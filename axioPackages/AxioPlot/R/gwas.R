plotGWASPValues.default <- function( position , chromosome , pValue , geneIDs , sort = TRUE , threshold = NULL , cutoff = NULL , suggestive = NULL , threshLabels = FALSE , geneLabels = FALSE , mainTitle = NULL , ... )
{
  if ( !is.null( cutoff ) )
  {
    coff <- pValue > -log10( cutoff ) & !is.na( pValue )
    position <- position[coff]
    chromosome <- chromosome[coff]
    pValue <- pValue[coff]
    geneIDs <- geneIDs[coff]
  }
  if ( sort )
  {
    ordered <- order(chromosome,position)
    position <- position[ordered]
    chromosome <- chromosome[ordered]
    pValue <- pValue[ordered]
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
  yLabel <- expression(-Log[10]~~group("(",P-Value,")"))
  if ( is.null( threshold ) )
  {
    yLim <- c( min(pValue,na.rm=TRUE)-0.1 , max(pValue)+0.1 )
  }
  else
  {
    yLim <- c( min(pValue,na.rm=TRUE)-0.1 , c( max(pValue)+0.1 , -log10( threshold )+0.1 )[which.max( c( max(pValue)+0.1 , -log10( threshold )+0.1 ) )] )
  }
  plot( NewPos , pValue , pch=20 , col = Color , main = mainTitle , xlab = xLabel , ylab = yLabel , ylim = yLim , axes = FALSE , ... )
  axis( 2 )
  axis( 1 , xTicks , labels = FALSE )
  axis( 1 , xNames , names( xNames ) , tick = FALSE )
  if ( !is.null( threshold ) )
  {
	sigGenes <- pValue > -log10( threshold )
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

plotGWASPValues.hexbin <- function( position , chromosome , pValue , geneIDs , sort = TRUE , threshold = NULL , cutoff = NULL , suggestive = NULL , threshLabels = FALSE , geneLabels = FALSE , mainTitle = NULL , ... )
{
  if ( !is.null( cutoff ) )
  {
    coff <- pValue > -log10( cutoff ) & !is.na( pValue )
    position <- position[coff]
    chromosome <- chromosome[coff]
    pValue <- pValue[coff]
    geneIDs <- geneIDs[coff]
  }
  if ( sort )
  {
    ordered <- order(chromosome,position)
    position <- position[ordered]
    chromosome <- chromosome[ordered]
    pValue <- pValue[ordered]
    geneIDs <- geneIDs[ordered]
  }
  NewPos <- cumsum( as.numeric( position ) )

  xTicks <- unique( chromosome )
  xTicks <- xTicks[ order( xTicks ) ]
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
  xTicks <- c( min( NewPos , na.rm = TRUE ) , xTicks )
  Color <- Color + 1

  if ( is.null( mainTitle ) )
  {
    mainTitle <- "GWAS Plot"
  }
  xLabel <- "Chromosome Map Position (MB)"
  yLabel <- expression(-Log[10]~~group("(",P-Value,")"))
  pushViewport( plotViewport( c( 5 , 4 , 2 , 2 ) ) )
  pushViewport( dataViewport( NewPos , pValue , name = "plotRegion" ) )
  grid.text( xLabel , y = unit( -3 , "lines" ) )
  grid.text( yLabel , x = unit( -3 , "lines" ) , rot = 90 )
  grid.text( mainTitle , y = unit( max( pValue ) + 0.5 , "native" ) , gp = gpar( fontface = "bold" , fontsize = 14 ) )
  yLim <- c( min(pValue,na.rm=TRUE)-0.1 , c( max(pValue)+0.1 , -log10( threshold )+0.1 )[which.max( c( max(pValue)+0.1 , -log10( threshold )+0.1 ) )] )
  hb1 <- hexbin( x = NewPos[Color == 1] , y = pValue[Color == 1] , ybnds = yLim )
  hb2 <- hexbin( x = NewPos[Color == 2] , y = pValue[Color == 2] , ybnds = yLim )
  grid.hexagons( hb1 , pen = "black" , style = "centroids" )
  grid.hexagons( hb2 , pen = "red" , style = "centroids" )
  grid.yaxis()
  grid.xaxis( at = xTicks , label = FALSE )
  grid.text( names( xNames ) , x = unit( xNames , "native" ) , y = unit( -1 , "lines" ) )
  if ( !is.null( threshold ) )
  {
    sigGenes <- pValue > -log10( threshold )
    if ( geneLabels )
    {
      if ( any( sigGenes ) )
      {
        grid.points( NewPos[sigGenes] , pValue[sigGenes] , gp = gpar( col = Color[sigGenes] ) , pch = 20 )
        grid.text( geneIDs[sigGenes] , x = unit( NewPos[sigGenes] , "native" ) , y = unit( pValue[sigGenes] , "native" ) , just = "left" )
      }
      else
      {
        warning( "No Genes Significant at p < " , threshold )
      }
    }
    grid.lines( x = unit( xTicks , "native" ) , y = unit( rep( -log10( threshold ) , length( xTicks ) ) , "native" ) , gp = gpar( col = "red" , lty = 3 ) )
    if ( threshLabels )
    {
      grid.text( paste( "(Genome-wide significant, p < " , threshold , ")" , sep = "" ) , x = unit( xTicks[2] , "native" ) , y = unit( -log10( threshold ) + 0.1 , "native" ) , just = "bottom" , gp = gpar( col = "blue" ) )
    }
  }
  if ( !is.null( suggestive ) )
  {
    grid.lines( x = unit( xTicks , "native" ) , y = unit( rep( -log10( suggestive ) , length( xTicks ) ) , "native" ) , gp = gpar( col = "red" , lty = 2 ) )
    if ( threshLabels )
    {
      grid.text( paste( "(Genome-wide suggestive, p < " , suggestive , ")" , sep = "" ) , x = unit( xTicks[2] , "native" ) , y = unit( -log10( suggestive ) + 0.1 , "native" ) , just = "bottom" , gp = gpar( col = "blue" ) )
    }
  }
}

plotGWASPValues <- function( position , chromosome , pValue , geneIDs , sort = TRUE , type = c( "points" , "hexbin" ) , threshold = NULL , cutoff = NULL , suggestive = NULL , threshLabels = FALSE , geneLabels = FALSE , mainTitle = NULL , ... )
{
  switch( type,
  "lines" = plotGWASValues.lines( position = position , chromosome = chromosome , value = pValue , geneIDs = geneIDs , sort = sort , threshold = threshold , cutoff = cutoff , suggestive = suggestive , threshLabels = threshLabels , geneLabels = geneLabels , mainTitle = mainTitle , ... ) ,
  "points" = plotGWASPValues.default( position = position , chromosome = chromosome , pValue = pValue , geneIDs = geneIDs , sort = sort , threshold = threshold , cutoff = cutoff , suggestive = suggestive , threshLabels = threshLabels , geneLabels = geneLabels , mainTitle = mainTitle , ... ) ,
  "hexbin" = plotGWASPValues.hexbin( position = position , chromosome = chromosome , pValue = pValue , geneIDs = geneIDs , sort = sort , threshold = threshold , cutoff = cutoff , suggestive = suggestive , threshLabels = threshLabels , geneLabels = geneLabels , mainTitle = mainTitle , ... ) ,
  "default" = plotGWASPValues.default( position = position , chromosome = chromosome , pValue = pValue , geneIDs = geneIDs , sort = sort , threshold = threshold , cutoff = cutoff , suggestive = suggestive , threshLabels = threshLabels , geneLabels = geneLabels , mainTitle = mainTitle , ... ) )
}
