# uniFilter.q

# Goal: choose informative genes that have 
#       high discriminatory power but low correlation with other selected genes

# Model expression levels of a given gene by a 
# one way analysis of variance model with heterogeneous variances
# rank features using test statistic of the hypothesis that all means are equal
# Use greedy strategy of removing features with correlation > corCut with current feature with highest BF statistic
# see Chen et. al., "Gene Selection for Multi-Class Prediction of Microarray Data"
# default is the Brown Forsythe statistic because Chen et. al. found it to perform the best

# testStatistic is the name of a function with arguments: feature, classes, nClasses, n
# see e.g. BrownForsytheFilter()
# other statistics considered by Chen et. al. are: 
# ANOVA F, Cochran, Welch test statistics

# NOTE: should standardize variables to have mean 0 and variance 1

#-------------------------------------------------------
# Arguments:
#-----------
#predMatrix, 
#                      y,  
#                      maxChooseN= ncol(predMatrix)/10, 
#                      corCut= .2, 
#                      testStatistic = BrownForsythe,
#                      standardize= T, 
#                      printIt= F
#-------------------------------------------------------
                      
uniFilter <- function(predMatrix, 
                      y,  
                      maxChooseN= ncol(predMatrix)/10, 
                      corCut= .2, 
                      testStatistic = BrownForsythe,
                      standardize= T, 
                      printIt= F)
{
 
	if(standardize) { # standardize variables to have mean 0 and variance 1
		predMatrix <- sweep(predMatrix, 2, colMeans(predMatrix), FUN="-")
		predMatrix <- sweep(predMatrix , 2, colStdevs(predMatrix), FUN="/")
	}

	featureIndices= 1:ncol(predMatrix)

	nClasses = table(y)
	n = length(y)
mem.tally.reset()
	# calculate Brown Forsythe statistic, used to rank features
	statistics <- apply(predMatrix, 2, testStatistic, classes= y, nClasses= nClasses, n= n)
mem.tally.report()
	
	chooseIndices = rep(NA, length= maxChooseN)
	chooseStatistics = rep(NA, length= maxChooseN)

	deleteIndices = numeric(0)  # initialize
	nSelected<- 0

	repeat {
		
		nSelected <- nSelected + 1
		cat(nSelected, "\t")

		if(length(deleteIndices) ) {
			featureIndices <- featureIndices[-deleteIndices]
			statistics <- statistics[-deleteIndices]
		}
		
		# stop when no more features to delete	
		if(length(featureIndices) <= 1) {
			if(length(featureIndices) == 1) {
				chooseIndices[nSelected] <- featureIndices
				chooseStatistics[nSelected] <- statistics
			}
			break
		}

		bestOneIndex = which(statistics == max(statistics))[1]  # if there are ties, take the first one

		bestOne <- chooseIndices[nSelected] <- featureIndices[bestOneIndex]
		chooseStatistics[nSelected] <- statistics[bestOneIndex]

		
		featureIndices= featureIndices[-bestOneIndex]
		statistics <- statistics[-bestOneIndex]

		# At this point, only calculate correlation with most recent top ranked gene
		# Need to investigate whether need to need to calculate correlation with all selected genes
		if(length(featureIndices) >= 0) {
			absCorrelationWithBest = abs(cor(predMatrix[,featureIndices], predMatrix[,bestOne]))
 
# The smaller the corCut, the more features deleted, and potentially the smaller the list
#			deleteIndices = which(absCorrelationWithBest > corCut)  
			deleteIndices = featureIndices[absCorrelationWithBest > corCut] 
			
		}

		if(printIt) {
			cat("Selected Gene: ", bestOne, "\n")

			cat("Deleted features: \t")
			cat(featureIndices[deleteIndices], "\n")

			cat("Remaining features: \t")		
			cat(featureIndices[featureIndices], "\n\n")
		}

		if(nSelected >= maxChooseN) break
mem.tally.report()

    }

	cat("\n\n")
	list(indices= chooseIndices[!is.na(chooseIndices)], statistics= chooseStatistics[!is.na(chooseStatistics)])
	
}

#BrownForsytheFilter(predMatrix= ALLData[,2:20], y= ALLData[,1], corCut= .6,  maxChooseN=10)
# BrownForsytheFilter(predMatrix= ALLData[,2:20], y= ALLData[,1],  maxChooseN=10)
# BrownForsytheFilter(predMatrix= ALLData[,2:20], y= ALLData[,1],  maxChooseN=10, printIt=T)
# BrownForsytheFilter(predMatrix= ALLData[,2:20], y= ALLData[,1],  corCut= .8, maxChooseN=19, printIt=T)


#dos.time(
#BFFiltered <- BrownForsytheFilter(predMatrix= ALLData[,2:ncol(ALLData)], y= ALLData[,1],  maxChooseN=1000)
#)



