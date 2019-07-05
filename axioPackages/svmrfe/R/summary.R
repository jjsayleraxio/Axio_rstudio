
summary.cvSvmRFE <- function( x , ... )
{
  cat( "Predicted Classes\n\n" )
  print( x$cvPredict )
  cat( "\nFeature List\n\n" )
  print( x$featureList )
  cat( "\nConfusion Matrix.\n Rows are true classes, columns are predicted classes.\n\n")
  print( round( x$cvConfusion , 3 ) )
}
