
# svmRFE.formula.q

#-------------------------------------------------------
# Arguments:
#-----------
# x   a formula object, with the response on the left of a ~ operator 
#     and the terms, separated by + operators, on the right. 
#     The response should be a factor. 

# data = a data frame in which to interpret the variables named in the formula

# subset  expression specifying which subset of observations should be used in 
#         the fit. This can be a logical vector (which is replicated to have 
#          length equal to the number of observations), a numeric vector 
#          indicating the observation numbers to be included, or a character 
#          vector of the observation names that should be included. 
#          All observations are included by default. 

# na.action   a function to filter missing data. This is applied to the model.frame 
#             after any subset argument has been applied. The default is na.fail, 
#             which returns an error if any missing values are found.  

# ...        additional arguments to svmRFE()
#-------------------------------------------------------

svmRFE.formula = function(x, data = NULL, subset, na.action = na.fail, kernel = "linear", ...)
{
	call <- match.call()
	if(!inherits(x, "formula"))
		stop("method is only for formula objects")
	m <- match.call(expand.dots = FALSE)
  if ( is.matrix( eval( m$data , sys.parent() ) ) )
  {
		m$data <- as.data.frame( data )
  }
  m$... <- NULL
	names( m )[2] <- "formula"
  m[[1]] <- as.name( "model.frame" )
  m <- eval( m , sys.parent() )
  return( m )
  Terms <- attr(m, "terms")
  print( Terms )
  print( names( m ) )
	attr(Terms, "intercept") <- 0
	x <- model.matrix(Terms, m)
	y <- model.extract(m, response)
	cat( x , "\n" )
	ret <- svmRFE(x, y, kernel = kernel , subset = subset, na.action = na.action, ...)
	
	ret$call <- call
	ret$terms <- Terms
	if(!is.null(attr(m, "na.action")))
		ret$na.action <- attr(m, "na.action")
	oldClass(ret) <- "svmRFE.formula"
	return(ret)
}

