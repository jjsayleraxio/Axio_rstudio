# entropyNextIndex.q

# adaptively calculates the index 'iUnsorted' in the function nextUnsortedIndex()
#               see file featureSetCalcs.q

# This function is not normally called by users

# To understand the arguments nIntervals, entropyLimit, meanLimit:
# see the pseudo code in Figure 10, page 19, of
#  "Entropy-based gene ranking without selection bias for the predictive 
#      classification of microarray data", 
#  Cesare Furlanello , Maria Serafini , Stefano Merler  and Giuseppe Jurman 
# BMC Bioinformatics 2003, 4:54     doi:10.1186/1471-2105-4-54
# The electronic version of this article is the complete one and can be 
#    found online at: http://www.biomedcentral.com/1471-2105/4/54


#-------------------------------------------------------
# Arguments:
#-----------
#index        'iUnsorted' in the function nextUnsortedIndex()
#weights       svm weights
#nFeatures     number of features

#nIntervals = floor( sqrt( nFeatures ) )  
#entropyLimit = 0.5 * logb( nIntervals , 2 )
#meanLimit = 0.2
                              
#-------------------------------------------------------
                              
entropyNextIndex <- function( index , weights , nFeatures = length( weights ) , nIntervals = floor( sqrt( nFeatures ) ) ,  
                              entropyLimit = 0.5 * logb( nIntervals , 2 ) , meanLimit = 0.2 )
{
  Zero <- sum( weights == 0 )#, na.rm = TRUE )
	if ( Zero )
	{
	  return( max( 1 , index - Zero ) )
	}
# linearly transform weights so they range between 0 and 1
	projectedWeights <- ( weights - min( weights ) ) / ifelse( max( weights ) == min( weights ) , 1 ,  max( weights ) - min( weights ) )
# split [0,1] into nIntervals intervals
	boundaries <- seq( 0 , 1 , 1 / nIntervals )
	p <- sapply( seq( 1 , nIntervals ) , function( i , bounds , y ) sum( bounds[i] <= y & y <= bounds[i+1] ) ,
               bounds = boundaries , y = projectedWeights ) / nFeatures
	p[p < .Machine$double.eps] <- .Machine$double.eps  # avoid logb(0,2)
	entropy <- -sum( p * logb( p , 2 ) , na.rm = TRUE )
	
	meanProjectedWeights <- mean( projectedWeights )
	if ( ( entropy > entropyLimit ) && ( meanProjectedWeights > meanLimit ) )
	{
	  omitLength <- sum( projectedWeights <= 1 / nIntervals )
	  return( index - omitLength )
	}
  else
  {
    projectedWeights[projectedWeights < .Machine$double.eps] <- .Machine$double.eps  # avoid logb(0,2)
    logProjWeights <- logb( projectedWeights , 2 )
    meanlogProjWeights <- mean( logProjWeights )
    A <- sum( logProjWeights <= meanlogProjWeights , na.rm = TRUE ) * 0.5
    beta <- 0
    while ( beta <= A )
    {
# Doesn't seem right, because meanlogProjWeights increases, therefore including more.  
# but want to include less
     	meanlogProjWeights = 0.5 * meanlogProjWeights
      beta <- sum( logProjWeights <= meanlogProjWeights , na.rm = TRUE )
    }
    return( max( 1 , index - beta ) )
  }
}
