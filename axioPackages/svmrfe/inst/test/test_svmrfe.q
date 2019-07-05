library("svmrfe")
module(ArrayAnalyzer)

design <- rbind(
c("A96ANIT0-1-LAD.CEL","A96ANIT0-1-LAD","0mg") ,
c("A96ANIT0-2-LA.CEL","A96ANIT0-2-LA","0mg") ,
c("A96ANIT0-3-LA.CEL","A96ANIT0-3-LA","0mg") ,
c("A96ANIT00.5-1-LA.CEL","A96ANIT00.5-1-LA","0.5mg") ,
c("A96ANIT00.5-2-LA.CEL","A96ANIT00.5-2-LA","0.5mg") ,
c("A96ANIT00.5-3-LB.CEL","A96ANIT00.5-3-LB","0.5mg") ,
c("A96ANIT01-1-LA.CEL","A96ANIT01-1-LA","1mg") ,
c("A96ANIT01-3-LA.CEL","A96ANIT01-3-LA","1mg") ,
c("A96ANIT05-1-LA.CEL","A96ANIT05-1-LA","5mg") ,
c("A96ANIT05-2-LA.CEL","A96ANIT05-2-LA","5mg") ,
c("A96ANIT05-3-LB.CEL","A96ANIT05-3-LB","5mg") ,
c("A96ANIT10-1-LD.CEL","A96ANIT10-1-LD","10mg") ,
c("A96ANIT10-2-LB.CEL","A96ANIT10-2-LB","10mg") ,
c("A96ANIT10-3-LA.CEL","A96ANIT10-3-LA","10mg") ,
c("B96ANIT0-1-L.CEL","B96ANIT0-1-L","0mg") ,
c("B96ANIT0-2-L.CEL","B96ANIT0-2-L","0mg") ,
c("B96ANIT0-3-L.CEL","B96ANIT0-3-L","0mg") ,
c("B96ANIT20-1-L.CEL","B96ANIT20-1-L","20mg") ,
c("B96ANIT20-2-L.CEL","B96ANIT20-2-L","20mg") ,
c("B96ANIT20-3-L.CEL","B96ANIT20-3-L","20mg") ,
c("B96ANIT20-4-L.CEL","B96ANIT20-4-L","20mg") ,
c("B96ANIT50-1-L.CEL","B96ANIT50-1-L","50mg") ,
c("B96ANIT50-2-L.CEL","B96ANIT50-2-L","50mg") ,
c("B96ANIT50-3-L.CEL","B96ANIT50-3-L","50mg") ,
c("B96ANIT50-4-L.CEL","B96ANIT50-4-L","50mg") ,
c("B96ANIT100-1-L.CEL","B96ANIT100-1-L","100mg") ,
c("B96ANIT100-2-L.CEL","B96ANIT100-2-L","100mg") ,
c("B96ANIT100-3-L.CEL","B96ANIT100-3-L","100mg") ,
c("B96ANIT100-4-L.CEL","B96ANIT100-4-L","100mg") )
design <- design[order(design[,1]),]
pdata <- new("phenoData",pData=data.frame(dosage=design[,3],rep=c(1:29)),varLabels=list( dosage="Dosage",rep="arbitrary number"))
exset <- ReadAffy(celfile.path = "/home/david/Insightful/Bioeffects/Data/AFRL-HEPB_ANIT_Liver 0.5-100/",phenoData = pdata)
retchar <- function(x)
{
  return( substring(x,70,nchar(x)) )
}
dimnames(exprs(exset))[[2]] <- apply(as.matrix(dimnames(exprs(exset))[[2]]),MARGIN=1,FUN=retchar)
nexset <- rma( exset )
boxplot(split(exprs(nexset),pData(nexset)[[2]]))
ff1 <- pOverA( 0.25 , log2(100) )
ff2 <- function(x) ( IQR(x) > 0.50 )
ff3 <- function(x) (median(2^x) > 300 )
subset <- genefilter(abs(exprs(nexset)),filterfun(ff1,ff2,ff3))
afsubs <- t(exprs(nexset)[subset,])
dosage <- pData(nexset)[,1]

#-----------------------------------------
# cross validation

# delete by half
tryafsubsCVHalf <- cvSvmRFE(afsubs, dosage)
wmf.graph("/home/david/Insightful/doc/cv_error.wmf", width=11,height=11)
plot(tryafsubsCVHalf)
dev.off()
rfe.barplotcv(tryafsubsCVHalf, 32)

library(annotate)
library(XML)
library(rae230acdf)
library(rae230aAnnoData)
genenames <- ((dimnames(afsubs)[[2]])[tryafsubsCVHalf$featureList])[1:32]
ratgenes <- pubmed(genenames)
absts <- getPMInfo(ratgenes)

# fun of svm.weight
# warray: table of weights for diff between class i,j

# fun of orderFeatures
# sumabsw
# entropy delete strategy
# DAH: does not complete, problem somewhere within single fold run
tryafsubsCVEntropy <- cvSvmRFE(afsubs, dosage, delete="entropy")

#TODO: FIX
plot(tryafsubsCVEntropy)
rfe.plotcv(tryafsubsCVEntropy)
rfe.barplotcv(tryafsubsCVEntropy,10)

# Test switchToDeleteOne = 10 with halving strategy
# DAH: problem in storing errorCV.  Stores double the number of errorCV values
tryafsubsSwitchHalf <- cvSvmRFE(afsubs, dosage, switchToDeleteOne = 10)

# Test switchToDeleteOne = 10 with halving strategy
tryafsubsEntropySwitchEntropy <- cvSvmRFE(afsubs, dosage, , delete="entropy", switchToDeleteOne = 10)

#make prediction 
#remeber to change ALLData features as final model feature set
tryafsubsSwitchHalf <- svmRFE(afsubs, dosage,returnModels = T)
# svm.predict <- predict(tryALLSwitchHalf$testedModels[[4]]$model,ALLData[,tryALLSwitchHalf$testedModels[[4]]$featureSet][,-2])
svm.predict <- predict(tryafsubsSwitchHalf$testedModels[[1]]$model,afsubs)
confusion.svm = table(svm.predict, dosage)
confusion.svm

tryafsubsSwitchHalf <- svmRFE(afsubs, dosage,returnModels = T)
#svm.predict <- predict(tryALLSwitchHalf$testedModels[[4]]$model,ALLData[,tryALLSwitchHalf$testedModels[[4]]$featureSet][,-2])
svm.predict <- predict(tryafsubsSwitchHalf$testedModels[[1]]$model,afsubs)
confusion.svm = table(svm.predict, dosage)
confusion.svm


svm.only <- svm(afsubs, dosage, kernel="linear")
svm.only.predict <- predict(svm.only,afsubs)
confusion.svm = table(svm.only.predict, dosage)
confusion.svm
(sum(confusion.svm) - sum(diag(confusion.svm)))/sum(confusion.svm)  #  0.5, i.e. 50% correct
sum(diag(confusion.svm))/sum(confusion.svm) # percent correct

svm.only <- svm(afsubs, dosage)
svm.only.predict <- predict(svm.only,afsubs)
confusion.svm = table(svm.only.predict, dosage)
confusion.svm
(sum(confusion.svm) - sum(diag(confusion.svm)))/sum(confusion.svm)  #  0.5, i.e. 50% correct
sum(diag(confusion.svm))/sum(confusion.svm) # percent correct

#
all.equal(tryafsubsSwitchHalf$testedModels[[1]]$model, svm.only)

library(svmrfe)

generateClassificationDataSet <- function( nTestFeatures = 128 , nTestSamples = 100 , nTestImportantFeatures = 10 , nTestCategories = 4 , errorMean = 100 , errorVariance = 12 , testLowerBound = 150 , testUpperBound = 200 )
{
  testCategorical <- sample( rep( paste( "Category" , c( 1:nTestCategories ) ) , nTestSamples / nTestCategories ) , nTestSamples )
  testCategorical <- as.factor( testCategorical )
  testData <- matrix( rnorm( nTestFeatures * nTestSamples , errorMean , errorVariance ) , ncol = nTestFeatures )
  testEffects <- unique( round( runif(nTestSamples,testLowerBound,testUpperBound) ) )
  for( i in 1:nTestImportantFeatures )
  {
    testEffectsTmp <- sample( testEffects , nTestCategories )
    for ( j in 1:nTestCategories )
    {
      testData[testCategorical == levels( testCategorical )[j],i] <- testData[testCategorical == levels( testCategorical )[j],i] + testEffectsTmp[j]
    }
  }
  return( list( classes = testCategorical , data = testData ) )
}

svmRfeTestCVEntropy <- cvSvmRFE( svmRfeTestData , svmRfeTestCategorical, delete = "entropy")
plot( svmRfeTestCVEntropy )
svmRfeTestCVEntropy

stratifiedSampleList.testSvmRfe <- c()
for(type in levels(svmRfeTestCategorical)) {
	stratifiedSampleList.testSvmRfe[[type]] = which(svmRfeTestCategorical==type)
}

nKeepFeatures = ncol(svmRfeTestData)
# decide how many variables to sample at each node
# recommended by Breiman in a tech. report
nRandomSplitVars <- logb(nKeepFeatures, base=2) + 1

# default for R implementation
nRandomSplitVars = sqrt(nKeepFeatures)
nRandomSplitVars = round(nRandomSplitVars)
nRandomSplitVars = 25

# set the number of trees in the forest
nTrees = 23
nTrees = 100

testSvmRfe.forest <- forest(
	svmRfeTestCategorical ~ .,
	data = data.frame(svmRfeTestData, svmRfeTestCategorical),
	classVote = T,
	nTrees = nTrees,
	boost = T,
	nRandomSplitVars = nRandomSplitVars,
	samplerControl = list(name="rfBootstrapSampling", indexList=stratifiedSampleList.testSvmRfe),
	cumulativeTraceOOB = TRUE,
	control= rfEachTreeControl(minsplit=3, minbucket=1))	

	confusionForest = table(tdosage, combined.forest$predictOOB)
confusionForest

varImportanceSplitImprove.RF <- rfVariableImportanceSplitImprove(testSvmRfe.forest)

rfOOBPrediction(testSvmRfe.forest, permuteVars=F)
z <- rfOOBPrediction(testSvmRfe.forest, permuteVars=T)
varImportancePermute.RF <- rfVariableImportancePermute(z)

testSvmRfe.Forest.best <- combinedData[,varImportanceSplitImprove.RF>0.018]
features.forests <- dimnames( combinedForest.best )[[2]]

rfMisclassification(testSvmRfe.forest)
rfMisclassification(z)
