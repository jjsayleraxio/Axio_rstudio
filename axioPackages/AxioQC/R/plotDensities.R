plotDensities <- function( data , date = "Date" , type = "quarterly" , level = "100ng" , pathRoot = "./" , xlab = "Concentration" , ylab = "Density" , main = "Density by Quarter, 100ng" , subset = NULL , showLimits = NULL , plotBootD = FALSE , cumulative = FALSE , PNG = FALSE )
{
  main <- paste( main , level , sep = ", " )
  densities <- data[["densities"]]
  boundaries <- data[["boundaries"]]
  bootPercentiles <- data[["bootPercentiles"]]
  bootDensities <- data[["bootDensities"]]
  if ( !is.null( subset ) )
  {
    boundaries <- boundaries[subset]
    if ( length( subset ) > 1 )
    {
      subset <- c( subset , "All" )
    }
    densities <- densities[subset]
    bootPercentiles <- bootPercentiles[subset]
    if( !is.null( bootDensities ) )
    {
      bootDensities <- bootDensities[subset]
    }
  }
  if ( cumulative )
  {
    aIndex <- which( names( densities ) == "All" )
    aIndex <- ifelse( length( aIndex ) == 0 , 1 , -aIndex )
    densities <- densities[aIndex]
    aIndex <- which( names( bootPercentiles ) == "All" )
    aIndex <- ifelse( length( aIndex ) == 0 , 1 , -aIndex )
    bootPercentiles <- bootPercentiles[aIndex]
    if( !is.null( bootDensities ) )
    {
      aIndex <- which( names( bootDensities ) == "All" )
      aIndex <- ifelse( length( aIndex ) == 0 , 1 , -aIndex )
      bootDensities <- bootDensities[aIndex]
    }
    fLevels <- orderDate( names( densities ) )
  }
  else
  {
    fLevels <- c( orderDate( names( densities )[-which( names( densities ) == "All" )] ) , "All" )
  }
  
  xlimit <- findXLimits( densities )
  ylimit <- findYLimits( densities )
  
  if( PNG )
  {
    png( paste( pathRoot , paste( ifelse( cumulative , "cumulative" , "sequential" ) , paste( ifelse( is.null( subset ) , "alldensities" , gsub( " " , "_" , subset ) ) , level , "_quarterly.png" , sep = "" ) , sep = "_" ) , sep = "/" ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  }
  tmp <- data.frame( do.call( rbind , lapply( densities , function( x ) return( cbind( x$x , x$y ) ) ) ) , factor( rep( names( densities ) , each = length( densities[[1]]$x ) ) , levels = fLevels ) )
  names( tmp ) <- c( "Conc" , "density" , "Time" )
  denPlot <- qplot( Conc , density , data = tmp , geom = "line" , colour = Time , xlim = xlimit , ylim = ylimit , main = ifelse( cumulative , paste( "Cumulative" , main ) , main ) , xlab = xlab , ylab = ylab ) + scale_colour_discrete( guide = guide_legend( title = NULL ) ) + theme( legend.position = "top" )
  btmp <- data.frame( row.names = names( bootPercentiles ) , do.call( rbind , lapply( bootPercentiles , function( x , y1 , y2 ) return( cbind( x[1] , x[1] , y1 , y2 ) ) , ylimit[1] , 0.2 * ylimit[2] ) ) , factor( names( densities ) , levels = fLevels ) )
  names( btmp ) <- c( "x" , "xend" , "y" , "yend" , "Time" )
  denPlot <- denPlot + geom_segment( aes( x = x , xend = xend , y = y , yend = yend , colour = Time ) , data = btmp )
  btmp <- data.frame( row.names = names( bootPercentiles ) , do.call( rbind , lapply( bootPercentiles , function( x , y1 , y2 ) return( cbind( x[2] , x[2] , y1 , y2 ) ) , ylimit[1] , 0.2 * ylimit[2] ) ) , factor( names( densities ) , levels = fLevels ) )
  names( btmp ) <- c( "x" , "xend" , "y" , "yend" , "Time" )
  denPlot <- denPlot + geom_segment( aes( x = x , xend = xend , y = y , yend = yend , colour = Time ) , data = btmp )
  btmp <- data.frame( row.names = names( bootPercentiles ) , do.call( rbind , lapply( bootPercentiles , function( x , y1 , y2 ) return( cbind( x[3] , x[3] , y1 , y2 ) ) , ylimit[1] , 0.1 * ylimit[2] ) ) , factor( names( densities ) , levels = fLevels ) )
  names( btmp ) <- c( "x" , "xend" , "y" , "yend" , "Time" )
  denPlot <- denPlot + geom_segment( aes( x = x , xend = xend , y = y , yend = yend , colour = Time ) , data = btmp )
  if ( plotBootD )
  {
    lines( bootDensities[[1]] , col = Colors[1] , lty = 2 )
  }
  for ( i in names( densities )[-1] )
  {
    if ( plotBootD )
    {
      lines( bootDensities[[i]] , col = Colors[i] , lty = 2 )
    }
  }
  legendText <- makeLegendText( bootPercentiles )
  ylims <- rev( seq( 1.5 * median( ylimit ) , 0.975 * ylimit[2] , length.out = length( legendText ) ) )
  ttmp <- data.frame( x = xlimit[1] , y = ylims , label = legendText )
  denPlot <- denPlot + annotate( "text" , label = ttmp$label , x = ttmp$x , y = ttmp$y , hjust = 0 , size = 32 / length( legendText ) )
  if ( !is.null( showLimits ) )
  {
    if ( length( showLimits ) > 1 )
    {
      btmp <- data.frame( row.names = names( bootPercentiles ) , do.call( rbind , lapply( bootPercentiles , function( x , y ) return( cbind( x[1] , y ) ) , 0.2 * ylimit[2] ) ) , factor( names( densities ) , levels = fLevels ) )
      btmp <- rbind( btmp , data.frame( row.names = names( bootPercentiles ) , do.call( rbind , lapply( bootPercentiles , function( x , y ) return( cbind( x[2] , y ) ) , 0.2 * ylimit[2] ) ) , factor( names( densities ) , levels = fLevels ) ) )
      btmp <- rbind( btmp , data.frame( row.names = names( bootPercentiles ) , do.call( rbind , lapply( bootPercentiles , function( x , y ) return( cbind( x[3] , y ) ) , 0.1 * ylimit[2] ) ) , factor( names( densities ) , levels = fLevels ) ) )
      names( btmp ) <- c( "x" , "y" , "Time" )
      denPlot <- denPlot + geom_text( aes( x = x , y = y , colour = Time , label = paste( round( btmp$x , 3 ) ) ) , data = btmp , hjust = 0 , angle = 90 , show_guide = FALSE )
    }
    else
    {
      text( bootPercentiles[[1]][1] , 0.2 * ylimit[2] , offset = 0.4 , col = Colors[1] , labels = paste( round( bootPercentiles[[1]][1] , 3 ) ) , srt = 90 , adj = c( 0 , 0.5 ) )
      text( bootPercentiles[[1]][2] , 0.2 * ylimit[2] , offset = 0.4 , col = Colors[1] , labels = paste( round( bootPercentiles[[1]][2] , 3 ) ) , srt = 90 , adj = c( 0 , 0.5 ) )
    }
  }
  if ( PNG )
  {
    print( denPlot )
    dev.off()
  }
  invisible( denPlot )
}
