crossValErrors <- function ( trueLabels , cvFit ) 
{
  nFeatureSets <- cvFit$nSplits
  delete <- cvFit$delete
  nUniqueLabels <- length( table( trueLabels ) )
  predictionErrors <- matrix( NA , ncol = length( nFeatureSets ) , nrow = cvFit$nFold )
  errorConfusion <- vector( length = length( nFeatureSets ) , "numeric" )
  if ( inherits( trueLabels , "numeric" ) )
  {
    individualErrors <- vector( length = length( nFeatureSets ) , "numeric" )
  }
  else
  {
    individualErrors <- matrix( NA , ncol = length( nFeatureSets ) ,
        nrow = nUniqueLabels , dimnames = list( dimnames( table( trueLabels ) )[[1]] , NULL ) )
  }
  for ( i in 1:cvFit$nFold )
  {
    ii <- cvFit$groups[[i]]
    tmp <- collateErrorsNonEntropy( trueLabels , nFeatureSets , nUniqueLabels , cvFit , errorConfusion , individualErrors , ii )
    predictionErrors[i,] <- tmp$errorConfusion
  }
  tmp <- collateErrorsNonEntropy( trueLabels , nFeatureSets , nUniqueLabels , cvFit , errorConfusion , individualErrors )
  return( list( errorSE = sqrt( apply( predictionErrors , 2 , var ) / cvFit$nFold ) ,
          errorCV = colMeans( predictionErrors ) ,
          individualErrors = tmp$individualErrors ,
          errorConfusion = tmp$errorConfusion ) )
}
