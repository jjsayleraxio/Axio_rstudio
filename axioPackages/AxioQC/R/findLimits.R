findXLimits <- function( x )
{
  xlimit <- matrix( unlist( lapply( x , function( x ) range( x[["x"]] ) ) ) , nrow = length( x ) , byrow = TRUE )
  xlimit <- c( min( xlimit[,1] ) ,  max( xlimit[,2] ) )
  return( xlimit + c( -0.0005 , 0.0005 ) * xlimit )
}

findYLimits <- function( y )
{
  ylimit <- matrix( unlist( lapply( y , function( x ) range( x[["y"]] ) ) ) , nrow = length( y ) , byrow = TRUE )
  ylimit <- c( min( ylimit[,1] ) ,  max( ylimit[,2] ) )
  return( ylimit + c( -0.05 , 0.35 ) * ylimit )
}
