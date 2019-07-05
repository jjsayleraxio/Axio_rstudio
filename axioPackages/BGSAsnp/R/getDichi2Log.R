get.dichi2.log <- function( x , df , scale )
{
  nu <- df / 2
  return( nu * log( nu ) - log( gamma( nu ) ) + nu * log( scale ) - ( nu + 1 ) * log( x ) - ( nu * scale / x ) )
}
