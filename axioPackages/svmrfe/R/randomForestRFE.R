randomForestRFE <- function( x , y , delete = "half" , nFold = 10 , switchToDeleteOne = 1 , minFeatures = 1 , ... )
{
  if ( !inherits( x , "matrix" ) )
  {
    trainingData <- as.matrix( x )
  }
  nFeatures <- ncol( x )
  nLabels <- length( y )
  nFeatureSets <- numFeatureSets( nFeatures , delete , switchToDeleteOne )
  nFold <- trunc( nFold )
  if ( nFold < 2 )
  {
    stop( "nFold should be greater than or equal to 2" )
  }
  else if ( nFold > nLabels )
  {
    stop( "nFold should be less than or equal to the number of observations" )
  }
  featureList <- matrix( 0 , nFeatures , nFeatureSets )
  Splits <- vector( length = nFeatureSets )
  
# data frame of factor predictions, one for each feature set considered
  cvPredict <- data.frame( matrix( y , nLabels , ifelse( delete == "entropy" , nFeatures , nFeatureSets ) ) )
  iUnsorted <- nFeatures
  i <- 1
  splitList <- vector( mode = "numeric" )
  pb <- txtProgressBar( 1 , nFeatureSets , 1 , style = 3 )
  featureSet <- featureIndices <- c( 1:iUnsorted )
  varImportance <- vector( mode = "list" , length = nFeatureSets )
  predError <- vector( length = nFeatureSets )
  errorSE <- vector( length = nFeatureSets )
  if ( is.numeric( y ) )
  {
    indError <- vector( length = nFeatureSets )
  }
  else
  {
    indError <- matrix( nrow = length( table( y ) ) , ncol = nFeatureSets , dimnames = list( dimnames( table( y ) )[[1]] , NULL ) )
  }
  tDelete <- delete
  while( iUnsorted >= 1 )
  {
    if ( iUnsorted <= switchToDeleteOne )
    {
      tdelete = "one"
    }
    tmp <- data.frame( Response = y , x[,featureSet] )
    predError[i] <- errorest( Response ~ . , data = tmp , model = randomForest , estimator = ifelse( is.numeric( y ) , "cv" , "632plus" ) )$error
    rfTmp <- randomForest( Response ~ . , data = tmp , ... )
    varImportance[[i]] <- importance( rfTmp )[,1]
    cvPredict[,i] <- rfTmp$predicted
    if ( is.numeric( y ) )
    {
      indError[i] <- predError[i]
    }
    else
    {
      indError[,i] <- rfTmp$confusion[,"class.error"]
    }
    if ( is.numeric( y ) )
    {
      errorSE[i] <- sqrt( mean( rfTmp$mse ) / length( rfTmp$mse ) )
    }
    else
    {
      errorSE[i] <- mean( apply( rfTmp$err.rate , 2 , sd ) )
    }
    tmp <- orderFeaturesRF( varImportance[[i]] , featureIndices[1:iUnsorted] , iUnsorted , tDelete , splitList , i , ... )
    featureIndices <- tmp$featureIndices
    if ( iUnsorted <= minFeatures )
    {
      featureList[seq( 1 , length( tmp$featureIndices ) ),i] <- tmp$featureIndices
      Splits[i] <- max( tmp$featureIndices )
      break
    }
    featureSet <- tmp$featureSet
    featureList[seq( 1 , length( tmp$featureIndices ) ),i] <- tmp$featureIndices
    splitList <- tmp$splitList
    Splits[i] <- max( tmp$featureIndices )
    i <- tmp$i
    iUnsorted <- tmp$iUnsorted
    setTxtProgressBar( pb , i )
  }
  if ( delete == "entropy" )
  {
    cvPredict <- cvPredict[,c( 1:i ),drop = FALSE]
    featureList <- featureList[,c( 1:i ),drop = FALSE]
    Splits <- Splits[c( 1:i )]
    predError <- predError[c( 1:i )]
    if ( is.numeric( y ) )
    {
      indError <- indError[c( 1:i )]
    }
    else
    {
      indError <- indError[,c( 1:i ),drop = FALSE]
    }
    errorSE <- errorSE[c( 1:i )]
  }
  setTxtProgressBar( pb , nFeatures )
  close( pb )
  nFeatures <- Splits
  bestIndex <- featureList[,max( which.min( predError ) )]
  bestIndex <- bestIndex[which( bestIndex > 0 )]
  bestModel <- randomForest( x = x[,bestIndex] , y = y , importance = TRUE , ...)
  if ( inherits( y , "numeric" ) )
  {
    confusionMatrix <- sum( ( cvPredict[,min( which.min( predError ) )] - y )^2 ) / ( length( y ) - 1 )
    errorCV <- apply( cvPredict , 2 , function( x , y ) sum( ( x - y )^2 ) / ( length( y ) - 1 ) , y )
    names( errorCV ) <- apply( featureList , 2 , function( x ) length( which( x != 0 ) ) )
    individualErrors <- rev( indError )
    nBoot <- NULL
  }
  else
  {
    confusionMatrix <- table( y , factor( cvPredict[,min( which.min( predError ) )] , levels = levels( y ) ) )
    errorCV <- predError
    individualErrors <- indError[,rev( seq( 1 , length( nFeatures ) ) )]
    nBoot <- 1000
  }
  cvSvmRFEOut <- list( errorCV = rev( errorCV ) ,
      errorSE = rev( errorSE ) ,
      individualErrors = indError ,
      nFeatures = rev( apply( featureList , 2 , function( x ) length( which( x != 0 ) ) ) ) ,
      featureList = featureList[,rev( seq( 1 , length( nFeatures ) ) ),drop = FALSE] ,
#      seed = seed ,
      cvPredict = cvPredict[,rev( seq( 1 , length( nFeatures ) ) ),drop = FALSE] ,
      cvConfusion = rev( predError ) ,
      bestModel = bestModel ,
      confusion = confusionMatrix ,
      numeric = inherits( y , "numeric" ) ,
      nBoot = nBoot ,
      nFold = nFold
  )
  oldClass( cvSvmRFEOut ) = "cvSvmRFE"
  cvSvmRFEOut
}
