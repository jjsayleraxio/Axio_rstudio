setClass( "powerObject" )

plot.powerObject <- function( x , y = NULL , SD = NULL , ... )
{
  legendLabels <- names( x@effects )
  legendColors <- seq( 1 , length( x@effects ) )
  powerPlot <- NULL
  samplePlot <- NULL
#  theme_set( powerTheme )
  
  if ( !any( c( is.null( x@power ) , is.null( x@sampleSize ) ) ) )
  {
    nContrasts <- length( x@effects )
    layoutMat <- matrix( c( rep( 1 , nContrasts ) , rep( 2 , nContrasts ) , seq( 3 , (nContrasts+2) ) ) , nrow = 3 , ncol = nContrasts , byrow = TRUE )
    layout( layoutMat )
  }
  else
  {
    nContrasts <- length( x@effects )
    layoutMat <- matrix( c( rep( 1 , nContrasts ) , seq( 2 , (nContrasts+1) ) ) , nrow = 2 , ncol = nContrasts , byrow = TRUE )
    layout( layoutMat )
  }
  if ( !is.null( x@power ) )
  {
    mainTitle <- paste( "Power to detect" , x@type , "genotype effect" )
    if ( x@stype == "Normal" )
    {
      subTitle <- paste( paste( "SD = ", ifelse( is.null( SD ) , x@SD , SD ) , sep = "" ) , paste( " , N = " , x@N , sep = "" ) , sep = "" )
    }
    else
    {
      subTitle <- paste( "N = " , x@N , sep = "" )
    }
    main <- paste( mainTitle , subTitle , sep = "\n" )
    xLabel <- "Minor Allele Frequency"
    yLabel <- "Power"
    tmp <- data.frame( maf = rep( x@MAF , ncol( x@power ) ) , power = as.numeric( matrix( x@power , ncol = 1 ) ) , group = factor( rep( legendLabels , each = nrow( x@power ) ) ) )
    powerPlot <- qplot( maf , power , data = tmp , group = group , col = group , geom = "line" , main = main , xlab = xLabel , ylab = yLabel ) + powerTheme()
    powerPlot <- powerPlot + scale_x_continuous( breaks = x@MAF , labels = x@MAF )
  }
  if ( !is.null( x@sampleSize ) )
  {
    mainTitle <- paste( "Sample size to detect" , x@type , "genotype effect" )
    if ( x@stype == "Normal" )
    {
      subTitle <- paste( paste( "SD = ", ifelse( is.null( SD ) , x@SD , SD ) , sep = "" ) , paste( " , Power = " , x@beta , sep = "" ) , sep = "" )
    }
    else
    {
      subTitle <- paste( "N = " , x@N , sep = "" )
    }
    main <- paste( mainTitle , subTitle , sep = "\n" )
    xLabel <- "Minor Allele Frequency"
    yLabel <- "Sample Size"
    tmp <- data.frame( maf = rep( x@MAF , ncol( x@sampleSize ) ) , sampleSize = as.numeric( matrix( x@sampleSize , ncol = 1 ) ) , group = factor( rep( legendLabels , each = nrow( x@sampleSize ) ) ) )
    samplePlot <- qplot( maf , sampleSize , data = tmp , group = group , col = group , geom = "line" , main = main , xlab = xLabel , ylab = yLabel ) + powerTheme()
    samplePlot <- samplePlot + scale_x_continuous( breaks = x@MAF , labels = x@MAF ) + guides( size = guide_legend() )
  }
  grid.newpage()
  print( powerPlot , vp = viewport( 1 , 0.35 , x = 0.5 , y = 0.8 ) )
  print( samplePlot , vp = viewport( 1 , 0.35 , x = 0.5 , y = 0.45 ) )
  xLabels <- c( "AA" , "Aa" , "aa" )
  xTicks <- seq( 1 , 3 )
  xLabel <- "Genotype"
  yLabel <- "Phenotype"
  vpLength <- 1 / length( x@effects )
  for ( i in seq( 1 , length( x@effects ) ) )
  {
    if ( inherits( x@effects[[i]] , "matrix" ) )
    {
      xMat <- t( x@effects[[i]] )
    }
    else
    {
      xMat <- x@effects[[i]]
    }
    if ( x@type == "main" )
    {
      tmp <- data.frame( Genotype = xTicks , Phenotype = as.numeric( xMat ) )
      pPlot <- qplot( Genotype , Phenotype , data = tmp , geom = "line" , main = names( x@effects )[i] ) + scale_x_continuous( breaks = xTicks , labels = xLabels )
    }
    else
    {
      if ( is.null( colnames( xMat ) ) )
      {
        eLev <- c( "Control" , paste( "Arm" , seq( 1 , ncol( xMat ) - 1 ) ) )
        eLev <- factor( rep( eLev , each = 3 ) , levels = c( "Control" , setdiff( eLev , "Control" ) ) )
      }
      else
      {
        eLev <- colnames( xMat )
        eLev <- factor( rep( eLev , each = 3 ) , levels = colnames( xMat ) )
      }
      tmp <- data.frame( Genotype = xTicks , Phenotype = as.numeric( xMat ) , Arm = eLev , Proportion = rep( x@proportion / 10 , each = 3 ) )
      pPlot <- qplot( Genotype , Phenotype , data = tmp , geom = "line" , main = names( x@effects )[i] , colour = Arm , size = Proportion , alpha = 0.5 ) + guides( size = "none" , alpha = "none" ) + scale_x_continuous( breaks = xTicks , labels = xLabels )
    }
    print( pPlot , vp = viewport( vpLength , 0.3 , x = ( i - 0.5 ) * vpLength , y = 0.15 ) )
  }
  grid.draw( outPlot <- grid.grab( warn = 1 , wrap = TRUE ) )
  invisible( outPlot )
}

if ( !isGeneric( "plot" ) )
{
  setGeneric( "plot" , function( x , y , ... ) standardGeneric( "plot" ) )
}

setMethod( "plot" , signature( x = "powerObject" , y = "missing" ) , plot.powerObject )

powerTheme <- function( base_size = 10 )
{
  theme_grey( base_size = base_size ) %+replace%
  theme(

          axis.line =            element_blank(),
          axis.text.x =          element_text(size = base_size * 0.6 , lineheight = 0.9, colour = "grey50", hjust = 1, angle = 90),
          axis.text.y =          element_text(size = base_size * 0.6, lineheight = 0.9, colour = "grey50", hjust = 1),
          axis.ticks =           element_line(colour = "grey50"),
          axis.title.x =         element_text(size = base_size),
          axis.title.y =         element_text(size = base_size, angle = 90),
          axis.ticks.length =    unit(0.15, "cm"),
          axis.ticks.margin =    unit(0.1, "cm"),
          
          legend.background =    element_rect( colour = NA ), 
          legend.key =           element_rect(fill = "grey95", colour = "white"),
          legend.key.size =      unit(0.5, "lines"),
          legend.text =          element_text(size = base_size * 0.6),
          legend.title =         element_blank() ,
          legend.position =      "bottom" ,
          legend.box =           "horizontal" ,
          legend.margin =        unit(0.05, "lines") ,
          
          panel.background =     element_rect(fill = "grey90", colour = NA), 
          panel.border =         element_blank(), 
          panel.grid.major =     element_line(colour = "white"),
          panel.grid.minor =     element_line(colour = "grey95", size = 0.25),
          panel.margin =         unit(0.05, "lines") ,
          
          strip.background =     element_rect(fill = "grey80", colour = NA), 
          strip.text.x =         element_text(size = base_size * 0.8),
          strip.text.y =         element_text(size = base_size * 0.8, angle = -90),
          
          plot.background =      element_rect(colour = NA),
          plot.title =           element_text(size = base_size),
          plot.margin =          unit(c(0.05, 0.25, 0.05, 0.05), "lines")
          )
}
