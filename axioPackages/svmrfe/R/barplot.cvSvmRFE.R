barplot.cvSvmRFE <- function( rfeCvRes , nbfeatures = NULL ,
                              uniqueFeatures = TRUE , col = "blue" , ... )
{
  if ( is.null( nbfeatures ) )
  {
  	nBestFeatures <- min( rfeCvRes$nFeatures[which( rfeCvRes$errorCV ==
  	                              min( rfeCvRes$errorCV ) )] )
  }
  nFold <- ncol( rfeCvRes$featureList )
  bestFeatures <- rfeCvRes$featureList[1:nBestFeatures,]
  freqBestFeatures <- sort(table(bestFeatures))
  if ( uniqueFeatures )
  {
    freqBestFeatures <- sort(table(bestFeatures)[1:nBestFeatures])
  }
  else
  {
    freqBestFeatures <- sort(table(bestFeatures))
  }
  barplot(freqBestFeatures, names = names(freqBestFeatures), horiz = TRUE ,
          col = col , ... )
  if ( uniqueFeatures )
  {
    title( main = paste("Frequency of the best", nBestFeatures,
           "features from each of the", nFold, "folds of Cross Validation") )
  }
  else
  {
    title( main = paste("Frequency of the best", nBestFeatures,
           "features from each of the", nFold, "folds of Cross Validation\n",
           "with a total of", length(freqBestFeatures), "features") )
  }
}
