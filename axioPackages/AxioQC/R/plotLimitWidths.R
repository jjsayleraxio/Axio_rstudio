plotLimitWidths <- function( data , pathRoot = "./" , level = "100ng" , type = "quarterly" , main = "Interval Width by Quarter" , xlab = "Quarter" , ylab = "Interval Width" , cumulative = FALSE , PNG = FALSE )
{
  limitWidth <- data.frame( width = unlist( lapply( data[["bootPercentiles"]] , function( x ) return( diff( x[1:2] ) ) ) ) , quarter = names( data[["bootPercentiles"]] ) , stringsAsFactors = FALSE )
  limitWidth <- limitWidth[-which( limitWidth$quarter == "All" ),]
  limitWidth$quarter <- factor( limitWidth$quarter , levels = limitWidth$quarter )
  limitWidthOut <- lm( limitWidth$width ~ as.numeric( limitWidth$quarter ) )
  ylimit <- range( limitWidth$width )
  ylimit <- ylimit + c( -0.05 , 0.1 ) * ylimit
  
  if ( PNG )
  {
    png( paste( pathRoot , paste( ifelse( cumulative , "cumulative_intervalwidth" , "intervalwidth" ) , level , type , "png" , sep = "." ) , sep = "/" ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  }
  limPlot <- qplot( quarter , width , data = limitWidth , geom = "point" , main = ifelse( cumulative , paste( "Cumulative" , paste( main , level , sep = ", " ) ) , paste( main , level , sep = ", " ) ) , xlab = xlab , ylab = ylab , ylim = ylimit , shape = I( 0 ) , solid = FALSE , size = I( 30 ) )
  tmp <- data.frame( x = seq( 1 , nrow( limitWidth ) ) , y = predict( limitWidthOut ) )
  limPlot <- limPlot + geom_line( aes( x = x , y = y ) , data = tmp , colour = "blue" )
  tmp <- data.frame( y = limitWidth$width , quarter = limitWidth$quarter , label = paste( "n =" , lapply( data[["densities"]] , function( x ) x$n ) )[seq( 1 , nrow( limitWidth ) )] )
  limPlot <- limPlot + geom_text( aes( label = label , x = quarter , y = y ) , data = tmp )
  limPlot <- limPlot + annotate( "text" , label = paste( "lm( width ~ quarter ), slope =" , paste( round( limitWidthOut$coefficients[2] , 2 ) , round( anova( limitWidthOut )[["Pr(>F)"]][1] , 2 ) , sep = ", p-value = " ) ) , x = 1 , y = ylimit[2] , hjust = 0 , colour = "blue" )
  if ( PNG )
  {
    dev.off()
  }
  invisible( limPlot )
}
