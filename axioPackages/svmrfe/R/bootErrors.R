bootErrors <- function( trainingData , trueLabels , cvFit , bootPlus = TRUE ) 
{
  nFeatureSets <- cvFit$nSplits
  delete <- cvFit$delete
  nUniqueLabels <- length( trueLabelTable <- table( trueLabels ) )
  looBootError <- rep( 0 , length = length( nFeatureSets ) )
  seBootError <- rep( 0 , length = length( nFeatureSets ) )
  if ( !is.numeric( trueLabels ) )
  {
    looIBootError <- matrix( 0 , nrow = nUniqueLabels , ncol = length( nFeatureSets ) )
    rownames( looIBootError ) <- names(  trueLabelTable )
    predClassProbs <- matrix( 0 , nrow = nUniqueLabels , ncol = length( nFeatureSets ) )
    rownames( predClassProbs ) <- names(  trueLabelTable )
    predIError <- matrix( 0 , nrow = nUniqueLabels , ncol = length( nFeatureSets ) )
    rownames( predIError ) <- names(  trueLabelTable )
  }
  else
  {
    bootPlus <- FALSE
    looIBootError <- rep( 0 , length = length( nFeatureSets ) )
    predIError <- rep( 0 , length = length( nFeatureSets ) )
  }
  for ( n in seq( 1 , length( nFeatureSets ) ) )
  {
    loobterr <- rep( 0 , length = length( trueLabels ) )
    featureIndx <- cvFit$featureList[which( cvFit$featureList[,n] > 0 ),n]
    for ( i in seq( 1 , length( trueLabels ) ) )
    {
      indx <- unlist( apply( cvFit$groups[[n]] , 2 , function( x , i ) which( !( i %in% x ) ) , i ) )
      for ( j in indx )
      {
        loobterr[i] <- loobterr[i] + ifelse(  is.numeric( trueLabels ) , ( predict( cvFit$model[[n]][[j]]$model , newdata = trainingData[i,featureIndx,drop = FALSE] ) - trueLabels[i] )^2 , as.numeric( predict( cvFit$model[[n]][[j]]$model , newdata = trainingData[i,featureIndx,drop = FALSE] ) != trueLabels[i] ) )
      }
      loobterr[i] <- loobterr[i] / length( indx )
      if ( !is.numeric( trueLabels ) )
      {
        looIBootError[as.character( trueLabels[i] ),n] <- looIBootError[as.character( trueLabels[i] ),n] + loobterr[i]
        predClassProbs[as.character( trueLabels[i] ),n] <-predClassProbs[as.character( trueLabels[i] ),n] + 1
      }
    }
    seBootError[n] <- sd( loobterr ) / sqrt( length( trueLabels ) )
    looBootError[n] <- mean( loobterr )
    if ( is.numeric( trueLabels ) )
    {
      looIBootError[n] <- looBootError[n]
    }
    else
    {
      looIBootError[,n] <- looIBootError[,n] / trueLabelTable
    }
  }
  predError <- unlist( lapply( seq( 1 , length( nFeatureSets ) ) , function( n , cvFit , trueLabels )
          {
            return( sum( unlist( lapply( seq( 1 , length( cvFit$model[[n]] ) ) , function( j , n , cvFit , trueLabels )
                            {
                              if ( is.numeric( trueLabels ) )
                              {
                                return( sum( ( fitted( cvFit$model[[n]][[j]]$model ) - trueLabels[cvFit$groups[[n]][,j]] )^2 ) / length( trueLabels ) )
                              }
                              else
                              {
                                return( sum( as.numeric( fitted( cvFit$model[[n]][[j]]$model ) != trueLabels[cvFit$groups[[n]][,j]] ) ) / length( trueLabels ) )
                              }
                            } , n , cvFit , trueLabels ) ) ) / length( cvFit$model[[n]] ) )
          } , cvFit , trueLabels ) )
  if ( is.numeric( trueLabels ) )
  {
    iPredError <- looBootError
  }
  else
  {
    predClassProbs <- predClassProbs / ( ncol( cvFit$groups[[1]] ) * length( trueLabels ) )
    iPredError <- matrix( unlist( lapply( seq( 1 , length( nFeatureSets ) ) , function( n , cvFit , trueLabels , trueLabelTable )
                {
                  return( rowSums( matrix( unlist( lapply( seq( 1 , length( cvFit$model[[n]] ) ) , function( j , n , cvFit , trueLabels , trueLabelTable )
                                      {
                                        s <- table( fitted( cvFit$model[[n]][[j]]$model ) , trueLabels[cvFit$groups[[n]][,j]] )
                                        diag( s ) <- 0
                                        return( apply( s , 1 , sum ) / trueLabelTable )
                                      } , n , cvFit , trueLabels , trueLabelTable ) ) , nrow = nUniqueLabels ) / length( cvFit$model[[n]] ) ) )
                } , cvFit , trueLabels , trueLabelTable ) ) , ncol = length( nFeatureSets ) )
  }
  sePredError <- unlist( lapply( seq( 1 , length( nFeatureSets ) ) , function( n , cvFit , trueLabels )
          {
            return( mean( unlist( lapply( seq( 1 , length( cvFit$model[[n]] ) ) , function( j , n , cvFit , trueLabels )
                            {
                              if ( is.numeric( trueLabels ) )
                              {
                                return( sd( ( fitted( cvFit$model[[n]][[j]]$model ) - trueLabels[cvFit$groups[[n]][,j]] )^2 ) / sqrt( length( trueLabels ) ) )
                              }
                              else
                              {
                                return( sd( as.numeric( fitted( cvFit$model[[n]][[j]]$model ) != trueLabels[cvFit$groups[[n]][,j]] ) ) / sqrt( length( trueLabels ) ) )
                              }
                            } , n , cvFit , trueLabels ) ) ) )
          } , cvFit , trueLabels ) )
  if ( bootPlus )
  {
    gammaHat <- colSums( matrix( trueLabelTable / sum( trueLabelTable ) , nrow = length( trueLabelTable ) , ncol = ncol( predClassProbs ) ) * ( 1 - predClassProbs ) )
    Rhat <- ( looBootError - predError ) / ( gammaHat - predError )
    wHat <- 0.632 / ( 1 - 0.368 * Rhat )
    return( list( errorSE = wHat * seBootError + ( 1 - wHat ) * sePredError ,
            errorCV = looBootError ,
            errorCVSE = seBootError ,
            individualErrors = wHat * looIBootError + ( 1 - wHat ) * iPredError ,
            errorBoot = predError ,
            errorSEBoot = sePredError ,
            errorConfusion = wHat * looBootError + ( 1 - wHat ) * predError ) )
  }
  else
  {
    return( list( errorSE = 0.632 * seBootError + 0.368 * sePredError ,
            errorCV = looBootError ,
            errorCVSE = seBootError ,
            individualErrors = 0.632 * looIBootError + 0.368 * iPredError ,
            errorBoot = predError ,
            errorSEBoot = sePredError ,
            errorConfusion = 0.632 * looBootError + 0.368 * predError ) )
  }
}
