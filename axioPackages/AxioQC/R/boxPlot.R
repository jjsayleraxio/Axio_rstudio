boxPlot <- function( data , x = "x" , y = "y" , time , xlab , ylab , main , timeName )
{
  ggplotCall <- paste( "ggplot( data , aes( " , x , " , " , y , " ) ) + geom_boxplot( fill = \"steelblue\" , colour = \"steelblue\" , outlier.colour = \"blue\" , alpha = 0.15 ) + labs( title = \"" , main , "\" ) + scale_y_continuous(  name = \"" , ylab , "\" ) + scale_x_discrete( name = \"" , xlab , "\" ) + geom_boxplot( aes( fill = " , time , " ) , alpha = 0.625 , outlier.color = \"blue\" ) + scale_fill_discrete( name = \"" , timeName , "\" )", sep = "" )
  return( eval( parse( text = ggplotCall ) ) )
}

boxPlotWithSmoothLine <- function( controlData , dnarna = "DNA" , main = "Predicted Concentration versus Date" , lcl10 , lcl100 , ucl10 , ucl100 , dateMin = NULL , dateMax = NULL , annotationList = NULL )
{
  lcl10Idx <- which( controlData[[paste( dnarna , "10" , sep = "" )]]$"10ng" > lcl10 )
  ucl10Idx <- which( controlData[[paste( dnarna , "10" , sep = "" )]]$"10ng" < ucl10 )
  idx10 <- intersect( lcl10Idx , ucl10Idx )
  lcl100Idx <- which( controlData[[paste( dnarna , "100" , sep = "" )]]$"100ng" > lcl100 )
  ucl100Idx <- which( controlData[[paste( dnarna , "100" , sep = "" )]]$"100ng" < ucl100 )
  idx100 <- intersect( lcl100Idx , ucl100Idx )
  controlDat <- rbind( data.frame( ExpectedConcentration = "10ng/ul" , PredictedConcentration = controlData[[paste( dnarna , "10" , sep = "" )]][idx10,"10ng"] , Date = controlData[[paste( dnarna , "10" , sep = "" )]][idx10,"QC_Date"] ) , data.frame( ExpectedConcentration = "100ng/ul" , PredictedConcentration = controlData[[paste( dnarna , "100" , sep = "" )]][idx100,"100ng"] , Date = controlData[[paste( dnarna , "100" , sep = "" )]][idx100,"QC_Date"] ) )

  if ( !is.null( dateMin ) )
  {
    controlDat <- controlDat[which( controlDat$Date > dateMin ),]
  }
  if ( !is.null( dateMax ) )
  {
    controlDat <- controlDat[which( controlDat$Date < dateMax ),]
  }
  
  lmOut <- lm( PredictedConcentration ~ Date + ExpectedConcentration , data = controlDat )
  lmPred <- stats::predict( lmOut , newdata = controlDat , se = TRUE )
  controlDat$ucl <- lmPred$fit + 1.96 * lmPred$se.fit
  controlDat$lcl <- lmPred$fit - 1.96 * lmPred$se.fit

  years <- as.Date( paste( unique( format( controlDat$Date , "%Y" ) ) , "01-01" , sep = "-" ) )
  years <- years[order( years )]
  years <- years[-1]

  controlBoxPlot <- ggplot( controlDat , aes( Date , PredictedConcentration ) ) + geom_boxplot( aes( group = Date ) , outlier.colour = "blue" ) + scale_x_date( breaks = "1 week" , labels = date_format( "%d-%b" ) ) + theme( axis.text.x = element_text( angle = 45 , hjust = 1 ) ) + labs( y = "Predicted Concentration" , x = "Date" ) + facet_grid( ExpectedConcentration ~ . , scales = "free_y" ) + ggtitle( main )
  controlBoxPlot <- controlBoxPlot + geom_vline( xintercept = as.numeric( years ) , colour = "blue" , linetype = "longdash" )
  if ( length( years ) > 0 )
  {
    for( i in seq( 1 , length( years ) ) )
    {
      controlBoxPlot <- controlBoxPlot + geom_text( data = NULL , x = as.numeric( years[i] ) , y = 0.95 * ucl10 , label = format( years[i] , "%Y" ) , size = 3 , hjust = 0 , show_guide = FALSE , colour = "blue" )
    }
  }
  controlBoxPlot <- controlBoxPlot + geom_smooth( aes( ymin = lcl , ymax = ucl ) , data = controlDat , alpha = 0.25 , method = "loess" )
  if ( !is.null( annotationList ) )
  {
    for ( i in seq( 1 , nrow( annotationList ) ) )
    {
      controlBoxPlot <- controlBoxPlot + geom_vline( xintercept = as.numeric( annotationList[i,1] ) , colour = "blue" , linetype = "longdash" ) + geom_text( data = NULL , x = as.numeric( annotationList[i,1] ) , y = runif( 1 , 0.50 , 0.975 ) * ucl100 , label = paste( annotationList[i,2] , format( annotationList[i,1] , "%d%b%Y" ) , sep = ", " ) , size = 3 , hjust = 0 , show_guide = FALSE , colour = "blue" )    }
  }
  return( controlBoxPlot )
}
