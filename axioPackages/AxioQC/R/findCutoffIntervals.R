findCutoffIntervals <- function( x , limits )
{
  return( quantile( x , limits ) )
}
