
print.cvSvmRFE <- function( x , ... )
{
  cat( "Predicted Classes\n\n" )
  print( x$cvPredict )
  cat( "\nFeature List\n\n" )
  print( x$featureList )
  if ( !is.null( x$nBoot ) )
  {
    if ( x$numeric )
    {
      cat("\nBootstrap Sum of Squared Error.\n\n")
      print( round( x$cvConfusion , 4 ) )
      minError <- min( which( x$cvConfusion == min( x$cvConfusion ) ) )
      cat("\nBootstrap Sum of Squared Error of best model.\n\n")
      print( round( x$cvConfusion[minError] , 4 ) )
      cat("\nSum of Squared Error from Best Model.\n\n")
    }
    else
    {
      cat("\nnBootstrap confusion errors.\n\n")
      print( round( x$cvConfusion , 4 ) )
      minError <- min( which( x$cvConfusion == min( x$cvConfusion ) ) )
      cat("\nBootstrap confusion error of best model.\n\n")
      print( round( x$cvConfusion[minError] , 4 ) )
      cat("\nBootstrap confusion error by class of best model.\n\n")
      print( round( x$individualError[,minError] , 4 ) )
      cat("\nConfusion Matrix from Best Model.\n Rows are true classes, columns are predicted classes.\n\n")
    }
  }
  else
  {
    if ( x$numeric )
    {
      cat("\nCross Validation Sum of Squared Error.\n\n")
      print( round( x$errorCV , 4 ) )
      minError <- min( which( x$errorCV == min( x$cvConfusion ) ) )
      cat("\nCross Validation Sum of Squared Error of best model.\n\n")
      print( round( x$errorCV[minError] , 4 ) )
      cat("\nSum of Squared Error from Best Model.\n\n")
    }
    else
    {
      cat("\nCross validation confusion errors.\n\n")
      print( round( x$errorCV , 4 ) )
      minError <- min( which( x$errorCV == min( x$errorCV ) ) )
      cat("\nCross validation confusion error of best model.\n\n")
      print( round( x$errorCV[minError] , 4 ) )
      cat("\nCross validation confusion error by class of best model.\n\n")
      print( round( x$individualError[,minError] , 4 ) )
      cat("\nConfusion Matrix from Best Model.\n Rows are true classes, columns are predicted classes.\n\n")
    }
  }
  print( x$confusion )
}
