plotFails <- function( fails , main = "Fails by Well" , rows = 8 , columns = 12 )
{
  x <- seq( 1 , columns )
  y <- seq( 1 , rows )
  
  circRads <- vector( mode = "list" , length = rows * columns )
  
  for ( i in seq( 1 , columns ) )
  {
    for ( j in seq( 1 , rows ) )
    {
      circRads[[( i - 1 ) * rows + j]] <- c( l = x[i] - 0.5 , r = x[i] + 0.5 , b = y[j] - 0.5 , t = y[j] + 0.5 , x = 0.5 , y = 0.5 )
    }
  }
  
  failLevels <- table( fails )
  failLevels <- data.frame( rowAlpha = substr( names( failLevels ) , 1 , 1 ) , column = as.numeric( substr( names( failLevels ) , 2 , nchar( names( failLevels ) ) ) ) , Value = as.integer( failLevels ) )
  failLevels$row = as.numeric( unlist( lapply( failLevels$rowAlpha , function( x ) which( rev( LETTERS[1:rows] ) %in% x ) ) ) )
  
  plotOut <- ggplot( expand.grid( x , y ) , aes( x = Var1 , y = Var2 ) ) + labs( x = "Column" , y = "Row" )
  plotOut <- plotOut + lapply( circRads , FUN =  function( x ) annotation_custom( circleGrob( gp = gpar( fill = "white" , color = "black" ) , x = x[5] , y = x[6] ) , xmin =  x[1] , xmax =  x[2] , ymin = x[3] , ymax = x[4] ) ) + coord_fixed( ratio = 1 ) + scale_y_continuous( breaks = c( 1:rows ) , labels = rev( LETTERS[1:rows] ) , limits = c( 0.875 , rows + 0.13 ) ) + scale_x_continuous( breaks = c( 1:columns ) )
  plotOut <- plotOut + annotate( "text" , failLevels$column , failLevels$row , label = failLevels$Value , colour = "red" )
  plotOut <- plotOut + ggtitle( main )
  return( plotOut )
}

plotHistFails <- function( fails , rows = 8 , columns = 12 )
{
  failLevels <- data.frame( row = substr( fails , 1 , 1 ) , column = substr( fails , 2 , nchar( fails ) ) )
  
  plotOutCols <- ggplot( failLevels , aes( x = column ) ) + labs( x = "Column" , y = "Count" ) + geom_histogram( aes( y = ..count.. ) , binwidth = 1 ) + scale_x_discrete( breaks = c( 1:columns ) , labels = c( 1:columns ) , limits = c( 1:columns ) , drop = FALSE )
  plotOutRows <- ggplot( failLevels , aes( x = row ) ) + labs( x = "Row" , y = "Count" ) + geom_histogram( aes( y = ..count.. ) , binwidth = 1 ) + scale_x_discrete( breaks = LETTERS[1:rows] , labels = LETTERS[1:rows] , limits = LETTERS[1:rows] , drop = FALSE )
  plotOutCols <- plotOutCols + ggtitle( "Fails by Column" )
  plotOutRows <- plotOutRows + ggtitle( "Fails by Row" )
  return( list( columns = plotOutCols , rows = plotOutRows ) )
}
