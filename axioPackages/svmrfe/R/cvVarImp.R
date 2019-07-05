
# cvVarImp.q

# Gives variable importance defined as the number of times that a feature
# is selected in the different folds of cross validation

# Bar plot of variable importance defined as the number of times that a feature
# is selected in the n folds of cross validation

# TODO: provide reference

#-------------------------------------------------------
# Arguments:
#-----------
# object 	has class "cvSvmRFE" such as returned by cvSvmRFE()
# nBest		number of best genes 
#-------------------------------------------------------

cvVarImp <- function (object, nBest)
{
    bestGenes <- object$featureList[1:nBest, ]
    freq.bestGenes <- sort(table(bestGenes))

# Horizontal plot better for labeling
     barplot(freq.bestGenes, names=names(freq.bestGenes), horiz = TRUE)
	 title(main = paste("Frequency of the best", nBest, "genes in the n folds of CV")) 
}