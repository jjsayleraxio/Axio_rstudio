
collateErrorsEntropy <- function( trueLabels , nFeatureSets , nUniqueLabels , cvFit , errorConfusion , individualErrors )
{
  trueLabelTable <- table( trueLabels )
  if ( !( inherits( trueLabels , "numeric" ) ) )
  {
    levelsFit <- as.character( levels( trueLabels ) )
  }
  for ( i in 1:max( nFeatureSets ) )
  {
    if ( inherits( trueLabels, "numeric" ) )
    {
      tmpConfusion <- sum( ( cvFit$cvPredict[,i] - trueLabels )^2 ) / ( length( y ) - 1 )
      errorConfusion <- errorConfusion + tmpConfusion
      individualErrors[i] <- tmpConfusion
    }
    else
    {
      tmpConfusion <- table( trueLabels , factor( cvFit$cvPredict[,i] , levels = levelsFit , labels = levelsFit ) )
      rowSumsS <- rowSums( tmpConfusion )
      rowSumsS[rowSumsS < 1] <- 1
      errorConfusion <- errorConfusion + ( tmpConfusion / rowSumsS )
      diag( tmpConfusion ) <- 0
      individualErrors[,i] <- apply( tmpConfusion, 1, sum ) / trueLabelTable
    }
  }
  errorConfusion <- errorConfusion / rowSums( errorConfusion )
  if ( inherits(trueLabels, "numeric" ) )
  {
    for ( i in 1:nUniqueLabels )
    {
      errorConfusion[i,] <- errorConfusion[i,] / trueLabelTable[i]
    }
  }
  return( list( errorConfusion = errorConfusion , individualErrors = individualErrors ) )
}

alignEntropySplits <- function( nFeatureSets , splitList )
{
  maxList <- which( nFeatureSets == max( nFeatureSets ) )
  if ( length( maxList ) > 1 )
  {
  	maxList <- maxList[1]
  }
  alignedSets <- matrix( nrow = length( nFeatureSets ) , ncol = nFeatureSets[maxList] )
  alignedSets[maxList,] <- splitList[[maxList]]
  setSequence <- c( 1:length( nFeatureSets ) )[-maxList]
  for ( i in setSequence )
  {
  	currentList <- splitList[[i]]
  	for ( j in 1:nFeatureSets[i] )
  	{
  	  diffSet <- abs( currentList[j] - alignedSets[maxList,] )
  	  foundCommonSplit <- TRUE
  	  while( foundCommonSplit )
  	  {
  	  	commonSplit <- min( which( diffSet == min( diffSet , na.rm = TRUE ) ) )
  	  	if ( tmpTest <- is.na( alignedSets[i,commonSplit] ) )
  	  	{
  	  	  alignedSets[i,commonSplit] <- currentList[j]
  	  	}
  	  	else
  	  	{
  	  	  diffSet[commonSplit] <- NA
  	  	}
  	  	foundCommonSplit <- !tmpTest
  	  }
  	}
  }
  return( alignedSets )
}

matrixFromEntropy <- function( groups , nSamps , alignedSets ,
                               cvPredict ,trueLabels)
{
  predictions <- as.data.frame( matrix( NA , nrow = nSamps ,
                                ncol = ncol( alignedSets ) ) )
  if ( inherits(trueLabels, "factor" ) )
  {
    levelsFit <- as.character( levels(trueLabels) )
    for ( i in 1:nrow(alignedSets) )
    {
      tmp <- as.character( unlist( cvPredict[[i]] ) )
      for ( j in seq( along = levelsFit ) )
      {
        tmp[(!is.na(match(tmp,j)))] <- levelsFit[j]
      }
      for ( j in seq( along = groups[[i]] ) )
      {
        predictions[groups[[i]][j],!(is.na(alignedSets[i,]))] <-
          as.character( tmp[[j]] )
      }
    }
  }
  else
  {
    for ( i in 1:nrow(alignedSets) )
    {
  	  predictions[groups[[i]],which( !is.na(alignedSets[i,]) )] <-
  	    unlist( cvPredict[[i]] )
    }
  }
  return( predictions )
}

errorCVEntropy <- function( predictionError , cvFit )
{
  errorCV <- colMeans( predictionError , na.rm = TRUE )
  minIndices <- vector( "list" , length = cvFit$nFold )
  maxMinIndices <- vector( length = cvFit$nFold )
  minIndices <- which( min( errorCV ) == errorCV, arr.ind=TRUE )
  maxMinIndices <- max( minIndices )
  return( list( errorCV = errorCV , maxMinIndices = maxMinIndices ) )
}

formatErrorEntropy <- function( trueLabels , nGroups , alignedSets ,
                                predictedClass )
{
  nFittedColumns <- sum(!is.na(alignedSets))
  nTotalColumns <- length( alignedSets )
  columnsNA <- which( is.na(alignedSets) )
  columnsNotNA <- which( !is.na(alignedSets) )
  yTmp <- matrix( NA, ncol=nTotalColumns, nrow=nGroups )
  yTmp[,columnsNotNA] <- matrix( trueLabels, ncol=nFittedColumns, nrow=nGroups )
  predictionErrors <- apply(yTmp != predictedClass, 2, sum) / nGroups
  predictionErrors[columnsNA] <- NA
  return( predictionErrors )
}
