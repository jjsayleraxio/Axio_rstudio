
#plot.cvSvmRFE.q
# plus, helper function errorBars()

# TODO: need to test and complete as indicated below

# plot method for object of class "cvSvmRFE"

#-------------------------------------------------------
# Arguments:
#-----------
# x  object of class "cvSvmRFE"
# type  which cv plots to plot
#   "total" cross validated error vs. number of genes
#   "individual"  Class-Specific CV Errors
#   "both"  both total and individual
# cex  graphical parameter

#-------------------------------------------------------

plot.cvSvmRFE <- function( x , type = c( "both" , "total" , "individual" ) , cex = 1 , key = TRUE ,
    tot.main = "Overall Cross Validation Error" , ind.main = "Class Specific Cross Validation Errors" ,
    xLab = "Number of Features" , yLab = "Percent CV Error" )
{
  if ( !is.null( x$nBoot ) )
  {
    tot.main <- "Overall 0.632 Bootstrap Error"
    ind.main <- "Class Specific 0.632 Bootstrap Errors"
    yLab <- "Percent 0.632 Bootstrap Error"
  }
  if ( x$numeric )
  {
    if ( !is.null( x$nBoot ) )
    {
      yLimits <- c( 0.975 * min( x$cvConfusion ) , 1.025 * max( x$cvConfusion )  )
    }
    else
    {
      yLimits <- c( 0.975 * min( x$errorCV ) , 1.025 * max( x$errorCV )  )
    }
    yLab <- "Mean Square Error"
    .INDIVIDUAL <- FALSE
    .TOTAL <- TRUE
  }
  else
  {
    yLimits <- c( -0.025 , 1.025 )
    type <- match.arg( type )
    if (type == "both")
    {
      .TOTAL <- TRUE
      .INDIVIDUAL <- TRUE
    }
    else if ( type == "total" )
    {
      .TOTAL <- TRUE
      .INDIVIDUAL <- FALSE
    }
    else if ( type == "individual" )
    {
      .TOTAL <- FALSE
      .INDIVIDUAL <- TRUE
    }
    else
    {
      stop("type must be either \"both\", \"total\", or \"individual\"")
    }
  }
  nc <- length( x$nFeatures )
  p <- vector( mode = "list" )
  if (.TOTAL)
  {
    if ( !is.null( x$nBoot ) )
    {
      tmp <- data.frame( errorCV = x$cvConfusion , genes = seq( 1 , nc ) , errorSE = x$errorSE )
      minError <- min( which( x$cvConfusion == min( x$cvConfusion ) ) )
    }
    else
    {
      tmp <- data.frame( errorCV = x$errorCV , genes = seq( 1 , nc ) , errorSE = x$errorSE )
      minError <- min( which( x$errorCV == min( x$errorCV ) ) )
    }
    limits <- aes( ymax = errorCV + errorSE , ymin = errorCV - errorSE )
    p[[1]] <- ggplot( tmp , aes( x = genes , y = errorCV ) ) + geom_line( colour = "blue" ) +
        scale_x_continuous( name = xLab , breaks = tmp$genes , labels = x$nFeatures ) +
        scale_y_continuous( name = yLab , limits = yLimits ) +
        geom_vline( xintercept = minError , colour = "yellow" ) +
        geom_hline( yintercept = ifelse( !is.null( x$nBoot ) , x$cvConfusion[minError][1] , x$errorCV[minError][1] ) , colour = "red" ) +
        labs( title = tot.main )
    if ( !x$numeric && is.null( x$nBoot ) )
    {
      p[[1]] <- p[[1]] + geom_errorbar( aes( ymax = errorCV + errorSE , ymin = errorCV - errorSE ) , width = 0.25 )
    }
    if ( length( x$nFeatures ) > 10 )
    {
      p[[1]] <- p[[1]] + theme( axis.text.x = element_text( angle = 90 , hjust = 1 ) )
    }
  }
  if (.INDIVIDUAL)
  {
    cats <- rownames( x$individualErrors )
    tmp <- data.frame( individualErrors = as.vector( t( x$individualErrors ) ) , genes = rep( seq( 1 , nc ) , length( cats ) ) , cat = rep( cats , each = ncol( x$individualErrors ) ) )
    p[[2]] <- ggplot( tmp , aes( x = genes , y = individualErrors , colour = cat ) ) + geom_line() +
        scale_x_continuous( name = xLab , breaks = seq( 1 , nc ) , labels = x$nFeatures ) +
        scale_y_continuous( name = yLab , limits = yLimits ) +
        labs( title = ind.main ) + theme( legend.position = "top" , legend.direction = "horizontal" , legend.title = element_blank() )
    if ( length( x$nFeatures ) > 10 )
    {
      p[[2]] <- p[[2]] + theme( axis.text.x = element_text( angle = 90 , hjust = 1 ) )
    }
  }
  if ( length( p ) == 1 )
  {
    print( p[[1]] )
  }
  else
  {
    grid.newpage()
    pushViewport( viewport( layout = grid.layout( 2 , 1 ) ) )
    for ( i in seq( 1 , length( p ) ) )
    {
      print( p[[i]] , vp = viewport( layout.pos.row = i , layout.pos.col = 1 ) )
    }
  }
  invisible( p )
}
