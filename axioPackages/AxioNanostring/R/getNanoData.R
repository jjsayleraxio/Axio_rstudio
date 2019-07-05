getNanoData <- function( fileCon )
{
  a <- readLines( fileCon )
  a <- c( "<doc>" , a , "</doc>" )
  b <- xmlParse( a , asText = TRUE )
  r <- xmlRoot( b )
  outPut <- list()
  outPut$header <- read.csv( file = textConnection( xmlValue( r[[1]] ) ) , check.names = FALSE )
  outPut$sampleAttributes <- read.csv( file = textConnection( xmlValue( r[[2]] ) ) , check.names = FALSE )
  outPut$laneAttributes <- read.csv( file = textConnection( xmlValue( r[[3]] ) ) , check.names = FALSE )
  outPut$codeSummary <- read.csv( file = textConnection( xmlValue( r[[4]] ) ) , check.names = FALSE )
  return( outPut )
}
