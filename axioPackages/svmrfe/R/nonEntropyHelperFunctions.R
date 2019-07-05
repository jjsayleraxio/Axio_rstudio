
collateErrorsNonEntropy <- function( y , nFeatureSets , nc , cvFit , errorConfusion , individualErrors , Group = NULL )
{
  if ( is.null( Group ) )
  {
    Group <- seq( 1 , length( y ) )
  }
  if ( !inherits( y , "numeric" ) )
  {
    trueLabelTable <- as.numeric( table( y[Group] ) )
    if ( any( trueLabelTable == 0 ) )
    {
      trueLabelTable[which( trueLabelTable == 0 )] <- 1
    }
  }
  for ( i in 1:length( nFeatureSets ) )
  {
    tmp <- errorValues( y , cvFit = cvFit$cvPredict[,i] , Group = Group , trueLabelTable = trueLabelTable )
    if ( inherits( y , "numeric" ) )
    {
      errorConfusion[i] <- tmp$errorConfusion
      individualErrors[i] <- tmp$individualErrors
    }
    else
    {
      errorConfusion[i] <- tmp$errorConfusion
      individualErrors[,i] <- tmp$individualErrors
    }
  }
  return( list( errorConfusion = errorConfusion , individualErrors = individualErrors ) )
}

errorValues <- function( y , cvFit = NULL , pModel = NULL , tData = NULL , Group , trueLabelTable = NULL )
{
  if ( is.null( pModel ) )
  {
    preds <- cvFit[Group]
  }
  else
  {
    preds <- predict( pModel , newdata = tData[Group,] )
  }
  if ( inherits( y , "numeric" ) )
  {
    s <- sum( ( preds - y[Group] )^2 ) / ( length( y[Group] ) - 1 )
    errorConfusion <- s
    individualErrors <- s
  }
  else
  {
    cvY <- factor( preds , levels = levels( y[Group] ) )
    s <- table( y[Group] , cvY )
    errorConfusion <- mean( 1 - ( diag( s ) / trueLabelTable ) )
    diag( s ) <- 0
    individualErrors <- apply( s , 1 , sum ) / trueLabelTable
  }
  return( list( errorConfusion = errorConfusion , individualErrors = individualErrors ) )
}

trainingError <- function( trainingData , trueLabels , cvFit , sIndx )
{
  tmp <- lapply( seq( 1 , length( cvFit ) ) , function( i , td , tl , cv , si ) return( errorValues( y = tl , pModel = cv[[i]] , tData = td , group = si[,i] ) ) )
  print( tmp )
  return( tmp )
}
