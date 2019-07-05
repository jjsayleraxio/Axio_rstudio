# BrownForsytheFilter.q

# Model expression levels of a given gene by a 
# one way analysis of variance model with heterogeneous variances
# test statistic of the hypothesis that all means are equal

# see Chen et. al., "Gene Selection for Multi-Class Prediction of Microarray Data"

#-------------------------------------------------------
# Arguments:
#-----------
# feature = one column of a prediction matrix
# classes = factor variable of classifications
# nClasses = number of classes
# nObs = number of observations ( == length(feature) == length(classes) )
#-------------------------------------------------------

BrownForsythe <- function( feature , classes , nClasses , nObs )
{
  centeredClassMeans <- tapply( feature , classes , mean ) - mean( feature )
  numerator <- sum( nClasses * centeredClassMeans * centeredClassMeans )
  denominator <- sum( ( 1 - ( nClasses / nObs ) ) * tapply( feature , classes , var) )
  return( numerator / denominator )
}
