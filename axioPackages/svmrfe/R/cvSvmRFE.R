# cvSvmRFE.q

# cvSvmRFE(), plus helper functions
# crossValRFE()
# crossValErrors()

#-------------------------------------------------------
# Arguments:
#-----------
# trainingData          prediction matrix

# trueLabels            response

# nFold       number of folds k in k-fold cross validation

# minFeatures    minimum number of features

# delete:	Delete strategy.  Specifies how many features to delete at each
#           iteration: none, half, one, square root of the remaining features.
#           "none" means just fit an svm without eliminating any variables.

# TODO: need to implement some of, and test, strategies none, sqrt, entropy.  
#      None should be trivial
# sqrt: need to test and account for switchToDeleteOne in numFeatureSets() 
#       see file featureSetCalcs.q

# switchToDeleteOne: 	iteration at which to switch to deleting one variable 
#                       at a time from the delete strategy specified by
#                       'delete'.  Often it makes sense to delete quickly in the
#                       first iteration(s) because many weights are close to
#                       zero.  However, eventually more caution is warranted.
#                       Note: using 'delete' = "entropy" is an adaptive way to
#                       achieve this.

# seed:  	random number seed.  Cross validation involves random sampling. 
#           Use the seed to reproduce results.
#-------------------------------------------------------

cvSvmRFE <- function( trainingData , trueLabels , nFold = 10 ,  minFeatures = 1 ,
    delete = c( "half" , "none" , "one" , "sqrt" , "entropy" ) , nBoot = NULL ,
    switchToDeleteOne = 0 , seed = NULL , kernel = "linear" , ... ) 
{
  delete <- match.arg( delete )
  nFeatures <- ncol( trainingData )

  nFeatureSets <- numFeatureSets( nFeatures , delete , switchToDeleteOne )

  if ( !is.null( seed ) )
  {
    set.seed( seed )
  }
  
# For delete == "one", not always possible to return all models.  
# Therefore, returnModels = F and predictData = left out data in call to svmRFE
  crossValRFEOut <- crossValRFE( trainingData = trainingData ,
      trueLabels = trueLabels , nFold = nFold , nBoot = nBoot ,
      nFeatureSets = nFeatureSets , minFeatures = minFeatures ,
      delete = delete , switchToDeleteOne = switchToDeleteOne ,
      kernel = kernel , ... )
  cat( "Starting error calculations..." )
  if ( is.null( nBoot ) )
  {
    errors <- crossValErrors( trueLabels , crossValRFEOut )
  }
  else
  {
    errors <- bootErrors( trainingData , trueLabels , crossValRFEOut )
  }
  cat( "Done!\n" )

  nFeatures <- numFeatureBySets( nFeatures, delete , switchToDeleteOne = switchToDeleteOne , cvFit = crossValRFEOut )

  bestIndex <- crossValRFEOut$featureList[,max( which.min( errors$errorConfusion ) )]
  bestIndex <- bestIndex[which( bestIndex > 0 )]
  bestModel <- svm( trainingData[,bestIndex] , trueLabels , kernel = kernel , ...)
  if ( inherits( trueLabels , "numeric" ) )
  {
    confusionMatrix <- sum( ( crossValRFEOut$cvPredict[,min( which.min( errors$errorCV ) )] - trueLabels )^2 ) / ( length( trueLabels ) - 1 )
    individualErrors <- rev( errors$individualErrors )
  }
  else
  {
    confusionMatrix <- table( trueLabels , factor( crossValRFEOut$cvPredict[,min( which.min( errors$errorCV ) )] , levels = levels( trueLabels ) ) )
    individualErrors <- errors$individualErrors[,rev( seq( 1 , length( nFeatures ) ) )]
  }

  cvSvmRFEOut <- list( errorCV = rev( errors$errorCV ) ,
      errorSE = rev( errors$errorSE ) ,
      individualErrors = individualErrors ,
      nFeatures = nFeatures ,
      featureList = crossValRFEOut$featureList[,rev( seq( 1 , length( nFeatures ) ) )] ,
      seed = seed ,
      cvPredict = crossValRFEOut$cvPredict[,rev( seq( 1 , length( nFeatures ) ) )] ,
      cvConfusion = rev( errors$errorConfusion ) ,
      bestModel = bestModel ,
      confusion = confusionMatrix ,
      numeric = inherits( trueLabels , "numeric" ) ,
      nBoot = nBoot ,
      nFold = nFold
  )
  oldClass( cvSvmRFEOut ) = "cvSvmRFE"
  return( cvSvmRFEOut )
}

crossValRFE <- function( trainingData , trueLabels , nFold = 10 , nFeatureSets , switchToDeleteOne = NULL ,
    kernel = "linear" , delete = "half" , minFeatures = 1 , nBoot = NULL , ... ) 
{
#  call <- match.call()
  if ( !inherits( trainingData , "matrix" ) )
  {
    trainingData <- as.matrix( trainingData )
  }
  nFeatures <- ncol( trainingData )
  nLabels <- length( trueLabels )
  if ( is.null( nBoot ) )
  {
    nFold <- trunc( nFold )
    if ( nFold < 2 )
    {
      stop( "nFold should be greater than or equal to 2" )
    }
    else if ( nFold > nLabels )
    {
      stop( "nFold should be less than or equal to the number of observations" )
    }
    else if ( nFold == nLabels )
    {
      groups <- 1:nLabels
      leaveOut <- 1
    }
    else
    {  # create nFold groups of indices
      leaveOut <- round( nLabels / nFold )
      groups <- createGroups( leaveOut , nLabels , nFold )
    }
  }
  else if ( !is.null( nBoot ) )
  {
    sIndxs <- vector( mode = "list" , length = nFeatureSets )
    modelFits <- vector( mode = "list" , length = nFeatureSets )
    bootPredict <- vector( mode = "list" , length = nBoot )
    testPredict <- vector( mode = "list" , length = nBoot )
  }
# data frame of factor predictions, one for each feature set considered
  cvPredict <- data.frame( matrix( trueLabels , nLabels , ifelse( delete == "entropy" , nFeatures , nFeatureSets ) ) )
  featureList <- matrix( 0 , nFeatures , nFeatureSets )
  Splits <- vector( length = nFeatureSets )

  iUnsorted <- nFeatures
  i <- 1
  splitList <- vector( mode = "numeric" )
  pb <- txtProgressBar( 1 , nFeatureSets , 1 , style = 3 )
  featureSet <- featureIndices <- c( 1:iUnsorted )
  while ( iUnsorted >= 1 )
  {
    rfeFit <- vector( mode = "list" , length = ifelse( is.null( nBoot ) , nFold , nBoot ) )
    if ( iUnsorted <= switchToDeleteOne )
    {
      delete = "one"
    }
    if ( is.null( nBoot ) )
    {
      for ( j in 1:nFold )
      {
        tpredData <- trainingData[groups[[j]],featureSet,drop = FALSE]
        dimnames( tpredData ) <- list( rownames( trainingData )[groups[[j]]] , names( trainingData ) )
        rfeFit[[j]] <- svmRFE( trainingData[-groups[[j]],featureSet,drop = FALSE] , trueLabels[-groups[[j]]] , minFeatures = minFeatures ,
            predictData = tpredData , kernel = kernel , delete = delete , CV = TRUE , ... )
        cvPredict[groups[[j]],i] <- rfeFit[[j]]$prediction
      }
    }
    else
    {
      sIndxMat <- matrix( NA , nLabels , nBoot )
      for ( j in seq( 1 , nBoot ) )
      {
        preds <- data.frame( matrix( trueLabels , nLabels , ifelse( delete == "entropy" , nFeatures , nFeatureSets ) ) )
        sIndx <- sample( seq( 1 , nLabels ) , replace = TRUE )
        nsIndx <- seq( 1 , nLabels )[which( !( seq( 1 , nLabels ) %in% sIndx ) )]
        sIndxMat[,j] <- sIndx
        tpredData <- trainingData[nsIndx,featureSet,drop = FALSE]
        colnames( tpredData ) <- colnames( trainingData )[featureSet]
        trData <- trainingData[sIndx,featureSet,drop = FALSE]
        rownames( trData ) <- seq( 1 , nLabels )
        rfeFit[[j]] <- svmRFE( trData , trueLabels[sIndx] , minFeatures = minFeatures ,
            predictData = tpredData , kernel = kernel , delete = delete , CV = TRUE , ... )
        cvPredict[nsIndx,i] <- rfeFit[[j]]$prediction
      }
      sIndxs[[i]] <- sIndxMat
      modelFits[[i]] <- rfeFit
    }
    tmp <- orderFeatures( lapply( rfeFit , function( x ) return( x$model ) ) , featureIndices[1:iUnsorted] , iUnsorted , delete ,
        returnModels = TRUE , kernel = kernel , CV = TRUE , splitList , i , ... )
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
    setTxtProgressBar( pb , ifelse( delete =="entropy" , nFeatureSets - length( featureIndices ) , i ) )
  }
  if ( delete == "entropy" )
  {
    cvPredict <- cvPredict[,c( 1:i ),drop = FALSE]
    featureList <- featureList[,c( 1:i ),drop = FALSE]
    Splits <- Splits[c( 1:i )]
    if ( !is.null( nBoot ) )
    {
      sIndxs <- sIndxs[c( 1:i )]
      modelFits <- modelFits[c( 1:i )]
    }
  }
  setTxtProgressBar( pb , nFeatures )
  close( pb )
  return( list( cvPredict = cvPredict , 
          nFold = nFold ,
          groups = eval( parse( text = ifelse( is.null( nBoot ) , "groups" , "sIndxs" ) ) ) ,   
          model = eval( parse( text = ifelse( is.null( nBoot ) , "NULL" , "modelFits" ) ) ) ,   
          featureList = featureList ,
          nSplits = Splits ,
          delete = delete ,
          splitList = splitList ) )
}
