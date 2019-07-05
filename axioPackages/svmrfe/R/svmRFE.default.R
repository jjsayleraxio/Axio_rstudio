# svmRFE.default.q

svmRFE <- function( x, ... )
{
  UseMethod( "svmRFE" )
}

#-------------------------------------------------------
# Arguments:
#-----------

#x		matrix of predictors

#y		factor response

#minFeatures   number of features in the smallest feature set
                      
#delete   Delete strategy.  Specifies how many features to delete at each
#         iteration:
#           none, half, one, square root of the remaining features.  "none"
#           means just fit an svm without eliminating any variables.

#switchToDeleteOne iteration at which to switch to deleting one variable 
#                       at a time from the delete strategy specified by
#                      'delete'.  Often it makes sense to delete quickly in the
#                       first iteration(s) because many weights are close to
#                       zero.  However, eventually more caution is warranted.
#                       Note: using 'delete' = "entropy" is an adaptive way to
#                       achieve this.
 
#returnModels       if TRUE, return all svm models

#predictData        data to predict.  Mostly used when called by cvSvmRFE() in 
#                   performing n-fold cross validation
                            
#-------------------------------------------------------

svmRFE.default <- function( x , y , minFeatures = 1 , delete = c( "half", "none", "one", "sqrt", "entropy" ) , 
						                 switchToDeleteOne = 0 , returnModels = FALSE , predictData = NULL ,
                             kernel = "linear" , CV = FALSE , scale = FALSE , ... ) 
{
  delete <- match.arg( delete )
  iUnsorted <- ncol( x )
  # iUnsorted is the last index of unsorted features.  
  # i.e., 1:iUnsorted are unsorted; iUnsorted+1 to dim(x)[2] are sorted
  # featureIndices starts out with all features at each iteration of rfe, update
  # iUnsorted according to the feature deletion strategy
  featureIndices <- 1:iUnsorted
  nFeatureSets <- numFeatureSets( iUnsorted , delete , switchToDeleteOne )
  if ( returnModels )
  {
    testedModels <- vector( "list" , length = nFeatureSets )
  }
  if ( !is.null( predictData ) )
  {
    if ( inherits( predictData , c( "named" , "vector" ) ) )
    {
      predictData <- matrix( predictData , nrow = 1 )
    }
    nrowPredict <- nrow( predictData )
    prediction <- data.frame( matrix( y[1:nrowPredict] , nrowPredict , nFeatureSets ) )  # use y to get factors
    returnModels2 <- TRUE  # need to specify returnModels = TRUE in orderFeatures()
    Predict <- TRUE
  }
  else
  {
    returnModels2 <- returnModels
    Predict <- FALSE
  }
  if ( CV )
  {
    model <- svm( x , y , kernel = kernel , scale = scale , ...)
    ans <- list( prediction = predict( model , predictData ) , model = model , splitList = NULL ,
        featureList = featureIndices , nFeatureSets = nFeatureSets )
  }
  else
  {
    i <- 1
    splitList <- vector( mode = "numeric" )
    if ( delete == "entropy" )
    {
      tiUnsorted <- 0
    }
    while ( iUnsorted > 1 )
    {
      if ( iUnsorted <= switchToDeleteOne )
      {
        delete = "one"
      }
      model <- svm( x[,featureIndices[1:iUnsorted],drop = FALSE] , y , kernel = kernel , ...)
      tmp <- orderFeatures( model , featureIndices , iUnsorted , delete ,
          returnModels = returnModels2 , kernel = kernel , CV = CV , splitList , i , ... )
      featureIndices <- tmp$featureIndices
      if ( returnModels )
      {
        testedModels[[i]] <- list( model = tmp$model, featureSet = tmp$featureSet )
      }
      if ( Predict )
      {
        if ( nrow( predictData ) > 1 )
        {
          prediction[,i] <- predict( tmp$model, predictData[,tmp$featureSet] )
        }
        else
        {
          dPredNames <- dimnames( predictData )
          tpredictData <- matrix( predictData[,tmp$featureSet] , 1 )
          dimnames( tpredictData ) <- list( dPredNames[[1]] , dPredNames[[2]][tmp$featureSet] )
          prediction[,i] <- predict( tmp$model, tpredictData )
        }
      }
      if ( delete == "entropy" )
      {
        if( iUnsorted == tiUnsorted )
        {
          break
        }
        tiUnsorted <- iUnsorted
      }
      else
      {
        if( iUnsorted <= minFeatures )
        {
          break
        }
      }
      splitList <- tmp$splitList
      i <- tmp$i
      iUnsorted <- tmp$iUnsorted
    }
# need to keep sequence of iUnsorted for delete= "entropy"
    ans <- list(featureList = featureIndices)
    if ( delete != "entropy" )
    {
      ans$nFeatureSets <- nFeatureSets
    }
    else
    {
      ans$nFeatureSets <- i - 1
    }
    if(returnModels)
    {
      if ( delete != "entropy" )
      {
        ans$testedModels <- testedModels
      }
      else
      {
        ans$testedModels <- testedModels[[1:(i-1)]]
      }
    }
    if(Predict)
    {
      if ( delete != "entropy" )
      {
        ans$prediction <- prediction
      }
      else
      {
        ans$prediction <- prediction[,1:(i-1)]
      }
    }
    ans$splitList <- splitList
  }
  oldClass( ans ) <- "svmRFE"
  return( ans )
}

# Order features and return next iUnsorted
orderFeatures <- function( model , featureIndices , iUnsorted , delete , returnModels = FALSE ,
    kernel = "linear" , CV = FALSE , splitList = NULL , i = NULL , ... )
{
  svmWeights <- svmWeight( model , CV )
# Order variables from best to worst
  featureOrder <- rev( order( svmWeights ) )
  featureIndices[1:iUnsorted] <- featureIndices[1:iUnsorted][featureOrder]
  if ( delete == "entropy" )
  {
    iUnsorted <- nextUnsortedIndex( index = iUnsorted , delete = delete ,
        weights = svmWeights , nFeatures = length( svmWeights ) )
  }
  else
  {
    iUnsorted <- nextUnsortedIndex( index = iUnsorted , delete = delete )
  }
  splitList[i] <- iUnsorted
  i <- i + 1
  ans <- list( iUnsorted = iUnsorted , featureIndices = featureIndices , 
      featureSet = featureIndices[1:iUnsorted] , i = i , splitList = splitList )
  if( returnModels )
  {
    ans$model <- model
  }
  return( ans )
}

# Order features and return next iUnsorted for Random Forest
orderFeaturesRF <- function( rfWeights , featureIndices , iUnsorted , delete , splitList = NULL , i = NULL , ... )
{
# Order variables from best to worst
  featureOrder <- rev( order( rfWeights ) )
  featureIndices[1:iUnsorted] <- featureIndices[1:iUnsorted][featureOrder]
  if ( delete == "entropy" )
  {
    iUnsorted <- nextUnsortedIndex( index = iUnsorted , delete = delete ,
        weights = rfWeights , nFeatures = length( rfWeights ) )
  }
  else
  {
    iUnsorted <- nextUnsortedIndex( index = iUnsorted , delete = delete )
  }
  splitList[i] <- iUnsorted
  i <- i + 1
  ans <- list( iUnsorted = iUnsorted , featureIndices = featureIndices , 
      featureSet = featureIndices[1:iUnsorted] , i = i , splitList = splitList )
  return( ans )
}

standardizeWeights <- function( modelWeights , CV = FALSE )
{
  if ( CV )
  {
    modelWeights[which( modelWeights == 0 , arr.ind = TRUE )] <- .Machine$double.eps
    if ( ncol( modelWeights ) < 2 )
    {
      modelWeights <- mean( ( modelWeights / sqrt( sum( modelWeights^2 ) ) )^2 )
    }
    else
    {
      modelWeights <- t( apply( modelWeights , 1 , function( x ) x / sqrt( sum( x^2 ) ) ) )
      modelWeights <- modelWeights * modelWeights
      meanWeights <- colMeans( modelWeights )
      sdWeights <- apply( modelWeights , 2 , sd )
      modelWeights <- meanWeights / sdWeights
      if ( any( is.infinite( modelWeights ) ) )
      {
        modelWeights[which( is.infinite( modelWeights ) )] <- .Machine$double.eps
      }
      if ( any( is.nan( modelWeights ) ) )
      {
        modelWeights[which( is.nan( modelWeights ) )] <- .Machine$double.eps
      }
    }
    return( modelWeights )
  }
  else
  {
    return( colMeans( modelWeights ) )
  }
}

svmWeight <- function ( model , CV = FALSE ) 
{
  if ( CV )
  {
    return( standardizeWeights( do.call( rbind , lapply( model , svmWeight , FALSE ) ) , CV ) )
  }
  else
  {
    return( t( model$coefs ) %*% model$SV )
  }
}
