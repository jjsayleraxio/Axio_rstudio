#' Create model from variables
#'
#' Create a model to use for the GWAS by running a linear model on phenotypes identified in indepvars without the snp component
#'
#' @param x A dataframe of phenotypes, samples by row, phenotypes by column
#' @param depVar The name of the dependent variable in x. Defaults to NULL. If null, inputs random normalized values in for the depvar
#' @param indepVars Phenotypes to include in the model
#' @param pCutOff The pvalue cutoff to determine significant phenotypes
#' @param exclude Phenotypes to exclude from the model. defaults to NULL
#' 
#' @return lm object
#' 
#' @examples 
#' testDF<-data.frame(var1=c(1,1,1,1,2,2,2,2),var2=c(2,2,2,2,1,1,1,2),var3=c(2,3,2,4,3,4,2,2),var4=c(-1,-5,-2,-3,-4,-5,-1,-3))
#'
#' newModel<-findModel(x = testDF,depVar = "var1",indepVars=c("var2","var3","var4"), pCutOff=0.05)
#'
#' formula(newModel)==formula("var1~var2")
#'
#' @export
#' 
findModel <- function( x , depVar = NULL , indepVars , pCutOff = 0.10 , exclude = NULL )
{
  if ( is.null( depVar ) )
  {
    cat( "\n" )
    x <- data.frame( depVar = rnorm( nrow( x ) ) , x )
    modelIndepVars <- model.frame( paste( "depVar" , paste( indepVars , collapse = " + " ) , sep = "~" ) , data = x )
    x.lm <- lm( modelIndepVars , data = x )
    
    xeigs <- svd( x.lm$qr$qr )$d
    cat( "\nCondition Number =" , sqrt( max( xeigs ) / min( xeigs ) ) , "\n" )
    if ( length( attr( terms( formula( x.lm ) ) , "term.labels" ) ) > 1 )
    {
      cat( "\nVariance Inflation Factors\n" )
      print( HH::vif( formula( x.lm ) , as.data.frame( unclass( x.lm$model ) ) ) )
    }
    
    return( x.lm )
  }
  if ( !is.null( exclude ) )
  {
    indepVars <- indepVars[-which( indepVars %in% exclude )]
  }
  SigIndepVars <- vector( length = length( indepVars ) )
  PVals <- vector( length = length( indepVars ) )
  names( PVals ) <- indepVars
  if ( inherits( x[,depVar] , "numeric" ) )
  {
    for ( i in seq( 1 , length( indepVars ) ) )
    {
      #fixed, was using "get", it was failing
      SigIndepVars[i] <- ( PVals[i] <- anova( lm( as.formula(paste0( depVar,"~", indepVars[i] )) , data = x ) )$"Pr(>F)"[1] ) < pCutOff
    }
    
    cat( "Individual Independent Variable Significance\n" )
    print( matrix( PVals , ncol = 1 , dimnames = list( indepVars , "P-Value" ) ) )
    
    modelIndepVars <- model.frame( paste( depVar , paste( indepVars[SigIndepVars] , collapse = " + " ) , sep = "~" ) , data = x )
    
    cat( "\n" )
    print( anova( x.lm <- lm( modelIndepVars , data = x ) ) )
    
    xeigs <- svd( x.lm$qr$qr )$d
    cat( "\nCondition Number =" , sqrt( max( xeigs ) / min( xeigs ) ) , "\n" )
    if ( length( attr( terms( formula( x.lm ) ) , "term.labels" ) ) > 1 )
    {
      cat( "\nVariance Inflation Factors\n" )
      print( HH::vif( formula( x.lm ) , as.data.frame( unclass( x.lm$model ) ) ) )
    }
    return( x.lm )
  }
  else
  {
    if ( inherits( x[,depVar] , "character" ) )
    {
      x[,depVar] <- factor( x[,depVar] )
    }
    if ( length( levels( x[,depVar] ) ) > 2 )
    {
      for ( i in seq( 1 , length( indepVars ) ) )
      {
        mldata <- mlogit.data( x , varying = NULL , choice = depVar , shape = "wide" )
        SigIndepVars[i] <- ( PVals[i] <- ( summary( mlogit( get( depVar ) ~ 1 | get( indepVars[i] ) , data = mldata ) )$lratio )$p.value[1] ) < pCutOff
      }
      
      matrix( PVals , ncol = 1 , dimnames = list( indepVars , "P-Value" ) )
      
      modelIndepVars <- as.formula( paste( depvar , paste( indepVars[SigIndepVars] , collapse = " + " ) , sep = " ~ " ) )
      modelIndepVars.lm <- as.formula( paste( depvar , paste( indepVars[SigIndepVars] , collapse = " + " ) , sep = " ~ " ) )
      depVaranova <- as.data.frame( matrix( nrow = sum( SigIndepVars ) , ncol = 2 ) )
      names( depVaranova ) <- c( "LR Statistic" , "Pr(Chi)" )
      rownames( depVaranova ) <- indepVars[SigIndepVars]
      
      for ( i in which( SigIndepVars ) )
      {
        reducedVars <- SigIndepVars
        reducedVars[i] <- FALSE
        naIndex <- which( !is.na( x[,indepVars[i]] ) )
        modelIndepVarsR <- as.formula( paste( depVar , paste( indepVars[reducedVars] , collapse = " + " ) , sep = " ~ " ) )
        respOut <- anova( multinom( modelIndepVarsR , data = x[naIndex,] ) , multinom( modelIndepVars , data = x[naIndex,] ) )
        depVaranova[indepVars[i],"LR Statistic"] <- respOut$"LR stat."[2]
        depVaranova[indepVars[i],"Pr(Chi)"] <- respOut$"Pr(Chi)"[2]
      }
      
      x.lm <- lm( modelIndepVars.lm , data = x )
      xeigs <- svd( x.lm$qr$qr )$d
      cat( "\nCondition Number =" , sqrt( max( xeigs ) / min( xeigs ) ) , "\n" )
      if ( length( attr( terms( formula( x.lm ) ) , "term.labels" ) ) > 1 )
      {
        cat( "\nVariance Inflation Factors\n" )
        print( HH::vif( modelIndepVars , data = x ) )
      }
      return( multinom( modelIndepVars , data = x ) )
    }
    else
    {
      for ( i in seq( 1 , length( indepVars ) ) )
      {
        SigIndepVars[i] <- ( PVals[i] <- anova( glm( get( depVar ) ~ get( indepVars[i] ) , data = x , family = "binomial" ) , test = "Chisq" )$"Pr(>Chi)"[2] ) < pCutOff
      }
      
      cat( "Individual Independent Variable Significance\n" )
      print( matrix( PVals , ncol = 1 , dimnames = list( indepVars , "P-Value" ) ) )
      
      modelIndepVars <- model.frame( paste( depVar , paste( indepVars[SigIndepVars] , collapse = " + " ) , sep = "~" ) , data = x )
      
      cat( "\n" )
      print( anova( x.lm <- glm( modelIndepVars , data = x , family = "binomial" ) , test = "Chisq" ) )
      
      xeigs <- svd( x.lm$qr$qr )$d
      cat( "\nCondition Number =" , sqrt( max( xeigs ) / min( xeigs ) ) , "\n" )
      if ( length( attr( terms( formula( x.lm ) ) , "term.labels" ) ) > 1 )
      {
        cat( "\nVariance Inflation Factors\n" )
        cnames <- names( x )
        x <- cbind( x , as.numeric( x[,depVar] ) - 1 )
        names( x ) <- c( cnames , paste( depVar , "Numeric" , sep = "." ) )
        modelIndepVars <- model.frame( paste( paste( depVar , "Numeric" , sep = "." ) , paste( indepVars[SigIndepVars] , collapse = " + " ) , sep = "~" ) , data = x )
        x.lm1 <- lm( modelIndepVars , data = x )
        print( HH::vif( formula( x.lm1 ) , as.data.frame( unclass( x.lm1$model ) ) ) )
      }
      return( x.lm )
    }
  }
}

#' find correlated Phenotypes
#'
#' @param x A dataframe of phenotypes, samples by row, phenotypes by column
#' @param critVal keep if phenotype in resulting correlation matric if it has a correlation with another phenotype greater than this value
#' 
#' @return returns correlation matrix. keeps all phenotypes with a correlation above critVal
#' 
#' @examples 
#' testDF<-data.frame(var1=c(1,1,1,1,2,2,2,2),var2=c(2,2,2,2,1,1,1,2),var3=c(2,3,2,4,3,4,2,2),var4=c(-1,-5,-2,-3,-4,-5,-1,-3))
#'
#' findCorrelated(x = testDF, critVal=.7)
findCorrelated <- function( x , critVal = 0.7 )
{
  corX <- cor( x )
  idx <- unique( unlist( apply( corX - diag( ncol( x ) ) , 2 , function( x ) which( abs( x ) > critVal ) ) ) )
  idx <- idx[order( idx )]
  return( corX[idx,idx] )
}
