
# predict.svmRFE

#-------------------------------------------------------
# Arguments:
#-----------
# object   object of class "cvSvmRFE"
# newdata  data frame containing the values at which predictions are required. 

# TODO: check that the predictors referred to in the right side of formula(object) 
#       are present by name in newdata. 

# TODO: If newdata is missing, return fitted values 

#-------------------------------------------------------

predict.svmRFE <- function ( object , newdata ) 
{
  nFeatureSets <- length(object$testedModels) - 1
  prediction <- data.frame(matrix(0, dim(newdata)[1], nFeatureSets))
  
  for (j in 1:nFeatureSets)
  {
    tnewdata <- newdata[, (object$testedModels[[j]])$featureSet]
    prediction[, j] <- predict((object$testedModels[[j]])$model, newdata[, (object$testedModels[[j]])$featureSet])
  }
  names( prediction ) <- paste("Model",1:nFeatureSets)
  cat("dim(prediction) = ", dim(prediction), "\n")
  return(prediction)
}

predict.cvSvmRFE <- function( object , newdata ) 
{
  cat("Here\n")
#    nFeatureSets <- length(object$testedModels) - 1
#    prediction <- data.frame(matrix(0, dim(newdata)[1], nFeatureSets))
#
#    for (j in 1:nFeatureSets)
#    {
#	  tnewdata <- newdata[, (object$testedModels[[j]])$featureSet]
#	  prediction[, j] <- predict((object$testedModels[[j]])$model, newdata[, (object$testedModels[[j]])$featureSet])
#    }
#	names( prediction ) <- paste("Model",1:nFeatureSets)
#    return(prediction)
}
