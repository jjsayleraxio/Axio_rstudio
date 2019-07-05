# Calculate SNP z value
#' @export
get.snp.z <- function( snpDat , y )
{
  freq <- rowMeans( snpDat / 2 )
  if ( inherits( y , "factor" ) )
  {
    yLevels <- levels( y )
  }
  else
  {
    yLevels <- names( table( y ) )
  }
  r0 <- rowSums( snpDat[,which( y == yLevels[2] )] == 0 , na.rm =  TRUE )
  s0 <- rowSums( snpDat[,which( y == yLevels[1] )] == 0 , na.rm =  TRUE )
  r1 <- rowSums( snpDat[,which( y == yLevels[2] )] == 1 , na.rm =  TRUE )
  s1 <- rowSums( snpDat[,which( y == yLevels[1] )] == 1 , na.rm =  TRUE )
  r2 <- rowSums( snpDat[,which( y == yLevels[2] )] == 2 , na.rm =  TRUE )
  s2 <- rowSums( snpDat[,which( y == yLevels[1] )] == 2 , na.rm =  TRUE )
  r <- r0 + r1 + r2
  s <- s0 + s1 + s2
  n0 <- r0 + s0
  n1 <- r1 + s1
  n2 <- r2 + s2
  n <- n0 + n1 + n2
  zStat <- ( ( s * r1 - r * s1 ) + 2 * ( s * r2 - s2 * r ) ) / sqrt( r * s / n ) / sqrt( n * ( n1 + 4 * n2 ) - ( n1 + 2 * n2 )^2 )
  return( zStat )
}
