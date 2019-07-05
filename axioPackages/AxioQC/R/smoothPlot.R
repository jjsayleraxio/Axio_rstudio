smoothPlot <- function( data , x = "x" , y = "y" , wQCDate , time , xlab , ylab , xBreakMajor , xBreakMinor , main , timeName )
{
  ggplotCall <- paste( "qplot( " , x , " , " , y , ", data = data , geom = \"boxplot\" , group = interaction( factor( " , wQCDate , " ) , ", time , " ) , fill = " , time , " , outlier.colour = \"blue\" ) + geom_smooth( aes( group = " , time , " , colour = " , time , " ) , method = \"loess\" , show_guide = FALSE ) + scale_y_continuous(  name = \"" , ylab , "\" ) + scale_x_date( breaks = \"" , xBreakMajor , "\" , minor_breaks = \"" , xBreakMinor , "\" , name = \"" , xlab , "\" , labels = date_format( \"%d-%b-%Y\" ) ) + labs( title = \"" , main , "\" ) + theme( axis.text.x = element_text( size = 7 , angle = -90 , vjust = 0.25 ) ) + scale_fill_discrete( name = \"" , timeName , "\" ) + scale_colour_hue( guide = FALSE )", sep = "" )
  return( eval( parse( text = ggplotCall ) ) )
}
