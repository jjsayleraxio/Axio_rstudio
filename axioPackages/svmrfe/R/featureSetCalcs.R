# featureSetCalcs.q
# helper functions that calculate aspects of feature sets

# number of Feature sets for a delete strategy

#-------------------------------------------------------
# Arguments:
#-----------
# nVars   number of variables

# delete  Delete strategy.  Specifies how many features to delete at each iteration:
#           none, half, one, square root of the remaining features.  "none" means
#           just fit an svm without eliminating any variables.

# switchToDeleteOne  iteration at which to switch to deleting one variable 
#                       at a time from the delete strategy specified by 'delete'.
#                       Often it makes sense to delete quickly in the first iteration(s)
#                       because many weights are close to zero.  However, eventually
#                       more caution is warranted.  Note: using 'delete' = "entropy"
#                       is an adaptive way to achieve this.
#-------------------------------------------------------

numFeatureSets <- function( nVars , delete , switchToDeleteOne = NULL )
{		
  return( length( numFeatureBySets( nVars , delete , switchToDeleteOne ) ) )
}

# Number of features in each set for a delete strategy
# used by cvSvmRFE
#TODO: check order; may need to be reversed
numFeatureBySets <- function( nVars , delete , switchToDeleteOne = NULL , cvFit = NULL )
{
  if ( delete == "half" )
  {
    nBySets <- nVars
    while( nVars > 1 )
    {
      nBySets <- c( nBySets , ( nVars <- ceiling( nVars / 2 ) ) )
      if ( !is.null( switchToDeleteOne ) )
      {
        if ( nBySets[length( nBySets )] < switchToDeleteOne )
        {
          nBySets <- c( nBySets[-length( nBySets )] , seq( ifelse( switchToDeleteOne %in% nBySets , switchToDeleteOne - 1 , switchToDeleteOne ) , 1 ) )
          break
        }
      }
    }
    nBySets <- rev( nBySets )
  }
  else if ( delete == "one" )
  {
    nBySets <- seq( 1 , nVars )
  }
  else if (delete == "sqrt")
  {
    nBySets <- nVars
    while( nVars > 1 )
    {
      remove <- floor( sqrt( nVars ) )
      nVars <-nVars - remove
      nBySets <- c( nBySets , nVars )
      if ( !is.null( switchToDeleteOne ) )
      {
        if ( nVars < switchToDeleteOne )
        {
          nBySets <- c( nBySets[-length( nBySets )] , seq( ifelse( switchToDeleteOne %in% nBySets , switchToDeleteOne - 1 , switchToDeleteOne ) , 1 ) )
          break
        }
      }
    }
    nBySets <- rev( nBySets )
  }
  else if ( delete == "entropy" )
  {
    if ( is.null( cvFit ) )
    {
      nBySets <- seq( 1 , nVars )
    }
    else
    {
      nBySets <- rev( unlist( apply( cvFit$featureList , 2 , function( x ) return( length( which( x > 0 ) ) ) ) ) )
    }
  }
  else
  {
    stop( "Delete strategy, " , delete , ", not yet implemented" )
  }
  return( nBySets )
}

# calculate next index which divides the unsorted part of the feature vector from the sorted part
# depends on delete strategy
# ... are any extra arguments (in addition to index) possibly needed for the different strategies
nextUnsortedIndex <- function( index , delete , switchToDeleteOne = NULL , ... )
{
  if ( !is.null( switchToDeleteOne ) )
  {
    if ( index <= switchToDeleteOne )
    {
      return( index - 1 )
    }
  }
  if ( delete == "half" )
  { 
    index <- max( 1 , ceiling( index / 2 ) )
  }
  if ( delete == "one" )
  {
    index <- index - 1
  }
  if ( delete == "sqrt" )
  {
    index <- max( 1 , index - floor( sqrt( index ) ) )
  }
  if ( delete == "entropy" )
  {
    index <- entropyNextIndex( index , ... )
  }
  if ( !is.null( switchToDeleteOne ) )
  {
    if ( index < switchToDeleteOne )
    {
      return( switchToDeleteOne )
    }
  }
  return( index ) 
}
