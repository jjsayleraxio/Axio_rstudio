gwasQQPlot <- function( pValues , main = "QQ Plot" , abLine = TRUE , logTransform = FALSE , ... )
{
  unifValues <- seq( 1 , length( pValues ) )
  unifValues <- ( unifValues - 0.5 ) / length( unifValues )
  pValues <- pValues[order(pValues)]
  if ( logTransform )
  {
    pValues <- -log10( pValues )
    unifValues <- -log10( unifValues )
    yLabel <- expression(-Log[10]~~group("(",P-Value,")"))
  }
  else
  {
    yLabel <- "P-Value"
  }
  plot( unifValues , pValues , xlab = "Uniform(0,1) Quantiles" , ylab = yLabel , main = main , ... )
  if ( abLine )
  {
    abline( c( 0 , 0 ) , c( 1 , 1 ) , col = "red" )
  }
}