# calculates Gene value
#' @export
getGeneZ <- function( snpZ , geneSnpList )
{
  nu0 <- 0
  sc0 <- 1
  geneZ <- unlist( lapply( geneSnpList , makeBetas , snpZ , nu0 , sc0 ) )
  return( geneZ )
}

makeBetas <- function( x , snpZ , nu0 , sc0 )
{
  beta <- snpZ[x]
  p <- length( beta )
  if ( p < 3 )
  {
    geneZ <- mean( beta , na.rm = TRUE)
  }
  else
  {
    V <- mean( beta^2 , na.rm = TRUE)
    geneZ <- sign( beta[which.max( abs( beta ) )] ) * log(1 + ( nu0 * sc0 + V * p ) / ( nu0 + p - 2 ) )
  }
  names( geneZ ) <- NULL
  return( geneZ )
}
