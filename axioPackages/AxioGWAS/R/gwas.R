#' Update model to include interaction terms or replace parts of a function
#'
#' Input an existing model, and will update the model based on your inputs
#'
#' @param x data.frame where the colnames are found in the Model. Usually the phenotypes plus the genotypes of a snp (required)
#' @param Model the original formula to update. Must be a formula or model.frame (required)
#' @param interactionTerm One of the colnames of x, to be added to the model, and interact with the genotype
#' @param replace Remove the specified term from the model
#'
#' @import stats
#'
#' @return a model.frame object
#'
#' @examples
#' testDF<-data.frame(var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),var4=c(4,4,4,4),Add=c(0,0,0,0))
#' Model<-formula("var1~var2+var3+var4")
#'
#' newModel<-updateModel(testDF,Model,interaction="var2",replace="var4")
#'
#' identical(formula(newModel),formula("var1~var2+var3+Add+var2:Add"))
#'
#' @export
updateModel <- function( x , Model , interactionTerm = NULL , replace = NULL )
{
  if(!inherits(Model,c("formula","data.frame"))){
    stop("Model must be either a formula or a model.frame")
  }
  newForm <- formula( Model )
  
  if ( !is.null( replace ) )
  {
    if ( inherits( replace , "matrix" ) || inherits( replace , "data.frame" ) )
    {
      for ( i in nrow( replace ) )
      {
        newForm <- update( newForm , paste( "~ . -" , paste( replace[i,] , collapse = "+" ) , sep = "" ) )
      }
    }
    else
    {
      newForm <- update( newForm , paste( "~ . -" , paste( replace , collapse = "+" ) , sep = "" ) )
    }
  }
  
  newForm <- update( newForm , paste( "~. +" , "Add" ) )
  
  if ( !is.null( interactionTerm ) )
  {
    newForm <- update( newForm , paste( "~. +" , paste( interactionTerm , "Add" , sep = ":" ) ) )
  }

  return( newForm )
}

#' Update model to include interaction terms or replace parts of a function for RR
#'
#' Input an existing model, and will update the model for RR based on your inputs
#'
#' @param x data.frame where the colnames are found in the Model. Usually the phenotypes plus a column named "Add" for the genotypes of a snp (required)
#' @param Model the original formula to update. Must be a formula or model.frame (required)
#' @param interactionTerm One of the colnames of x, to be added to the model, and interact with the genotype
#' @param replace Remove the specified term from the model
#'
#' @import stats
#'
#' @return a list, with a model.frame object and a data.frame
#'
#' @examples
#' testDF<-data.frame(var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),var4=c(4,4,4,4),Add=c(0,0,0,0))
#' Model<-formula("var1~var2+var3+var4")
#'
#' newModel<-updateModelRR(testDF,Model,interactionTerm="var2",replace="var4")
#'
#' identical(formula(newModel$newForm),formula("var1~var2+var3+Add+var2:Add"))
updateModelRR <- function( x , Model , interactionTerm = NULL , replace = NULL )
{
  newForm <- formula( Model )
  if ( inherits( Model , "formula" ) )
  {
    Model <- model.frame( Model , x )
  }
  respVar <- names( attr( terms( Model ) ,"dataClasses" ) )[attr( terms( Model ) , "response" )]
  eval( parse( text = paste( paste( "x" , paste( respVar , "Num" , sep = "" ) , sep = "$" ) , "<-" , paste( "as.numeric( x$" , respVar , " ) - 1" , sep = "" ) ) ) )
  respVar <- paste( respVar , "Num" , sep = "" )
  if ( !is.null( replace ) )
  {
    if ( inherits( replace , "matrix" ) || inherits( replace , "data.frame" ) )
    {
      for ( i in nrow( replace ) )
      {
        newForm <- update( newForm , paste( "~ . -" , paste( replace[i,] , collapse = "+" ) , sep = "" ) )
      }
    }
    else
    {
      newForm <- update( newForm , paste( "~ . -" , paste( replace , collapse = "+" ) , sep = "" ) )
    }
  }
  if ( !is.null( interactionTerm ) )
  {
    newForm <- paste( respVar , paste( c( attr( terms( newForm ) , "term.labels" ) , paste( c( "Add" , paste( interactionTerm , "Add" , sep = ":" ) ) , collapse = " + " ) ) , collapse = " + " ) , sep = " ~ " )
  }
  return( list( newForm = newForm , x = x ) )
}
#' Update model to include Add as variable of a function
#'
#' Input an existing model, and will update the model with Add as a variable for the snp, strips interaction
#'
#' @param x data.frame where the colnames are found in the Model. Usually the phenotypes plus a column named "Add" for the genotypes of a snp (required)
#' @param Model the original formula to update. Must be a formula or model.frame (required)
#'
#' @return a model.frame object
#'
#' @examples
#' testDF<-data.frame(var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),var4=c(4,4,4,4),Add=c(0,0,0,0))
#' Model<-formula("var1~var2+var3+var4")
#'
#' newModel<-reduceModel(testDF,Model)
#'
#' identical(formula(newModel),formula("var1~var2+var3+Add"))
reducedModel <- function( x , Model )
{
  newForm <- formula( Model )
  if ( inherits( Model , "formula" ) )
  {
    Model <- model.frame( Model , x )
  }
  respVar <- names( attr( terms( Model ) ,"dataClasses" ) )[attr( terms( Model ) , "response" )]
  newForm <- paste( respVar , paste( c( attr( terms( newForm ) , "term.labels" ) , "Add" ) , collapse = " + " ) , sep = " ~ " )
  return( newForm )
}

#' Update model to include and 2 genotypes, interaction terms or replace parts of a function
#'
#' Input an existing model, and will update the model based on your inputs
#'
#' @param x data.frame where the colnames are found in the Model. Usually the phenotypes plus the genotypes of a snp (required)
#' @param Model the original formula to update. Must be a formula or model.frame (required)
#' @param interactionTerm One of the colnames of x, to be added to the model, and interact with the genotype
#' @param replace Remove the specified term from the model
#'
#' @import stats
#'
#' @return a model.frame object
#'
#' @examples
#' testDF<-data.frame(var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),var4=c(4,4,4,4),Add=c(0,0,0,0))
#' Model<-formula("var1~var2+var3+var4")
#'
#' newModel<-updateModel(testDF,Model,interaction="var2",replace="var4")
#'
#' identical(formula(newModel),formula("var1~var2+var3+Add+var2:Add"))
#'
#' @export
updateModel_2g <- function( x , Model , interactionTerm = NULL , replace = NULL )
{
  if(!inherits(Model,c("formula","data.frame"))){
    stop("Model must be either a formula or a model.frame")
  }
  newForm <- formula( Model )
  
  if ( !is.null( replace ) )
  {
    if ( inherits( replace , "matrix" ) || inherits( replace , "data.frame" ) )
    {
      for ( i in nrow( replace ) )
      {
        newForm <- update( newForm , paste( "~ . -" , paste( replace[i,] , collapse = "+" ) , sep = "" ) )
      }
    }
    else
    {
      newForm <- update( newForm , paste( "~ . -" , paste( replace , collapse = "+" ) , sep = "" ) )
    }
  }
  
  newForm <- update( newForm , paste( "~. +" , paste( c( "Add" , "Add2" , "Add2:Add" ) , collapse= " + " ) ) )
  
  if ( !is.null( interactionTerm ) )
  {
    newForm <- update( newForm , paste( "~. +" , paste( interactionTerm , "Add" , sep = ":" ) ) )
  }
  
  return( newForm )
}

#' Generate PGx data
#'
#' Input an Phenotype DF, snpdata, chromosome and gender var and output a single dataframe
#'
#' @param x data.frame of phenotypes
#' @param g vector of genotypes, In the format "#/#",
#' @param chromosome the chromosome the genotype is from
#' @param genderVar the colname of the gender variable in x. Make sure the values are either "Male" or "Female"
#'
#' @return a dataframe combining phenotype and genotype
#'
#' @examples
#' testDF<-data.frame(var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),gender=c("Male","Female","Male","Female"))
#' genotype<-c("A/A","A/B","A/A","B/B")
#'
#' pgx<-makePGxData(testDF,genotype,1,"gender")
#'
#' identical(pgx,data.frame(testDF,Add=c(0,1,0,2)))
#'
#' @family makePGxData
#' @export
makePGxData <- function( x , g , chromosome , genderVar,... ) {
  UseMethod("makePGxData",g)
}

#' @family makePGxData
#' @export
makePGxData.default<-function( x , g , chromosome , genderVar ){
  stop(paste("No method makePGxData for class:",class(g)[1])) 
}

#' @family makePGxData
#' @export
makePGxData.numeric<-function( x , g , chromosome , genderVar ){
  if ( min( g ) > 0 && max( g ) > 2 ){
      g <- g - min(g) # scale to 0
  }
  return( data.frame( x , Add = g ) )
}

#' @family makePGxData
#' @export
makePGxData.integer<-makePGxData.numeric

#' @family makePGxData
#' @export
makePGxData.character<-function( x , g , chromosome , genderVar ){
  Add <- makeAdditive( g , chromosome , ifelse(!is.null(genderVar),list(as.character( x[,genderVar] ) ),list(NULL))[[1]])
  return( data.frame( x , Add ) )
}

#' @family makePGxData
#' @export
makePGxData.data.frame<-function( x , g , chromosome , genderVar ){
  g <- unlist( c( g ) )
  return( makePGxData( x , g , chromosome , genderVar ) )
}


#' Generate LME PGx data
#'
#' Input an Phenotype DF, snpdata, chromosome and gender var and output a single dataframe
#'
#' @param x data.frame of phenotypes
#' @param g vector of genotypes, In the format "#/#",
#' @param chromosome the chromosome the genotype is from
#' @param genderVar the colname of the gender variable in x. Make sure the values are either "Male" or "Female"
#'
#' @return a dataframe combining phenotype and genotype
#'
#' @examples
#' testDF<-data.frame(ID=paste0("sub",1:4),var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),gender=c("Male","Female","Male","Female"))
#' genotype<-c(sub1="A/A",sub2="A/B",sub3="A/A",sub4="B/B")
#'
#' pgx<-makeLMEPGxData(testDF,genotype,1,"gender","ID")
#'
#' identical(pgx,data.frame(testDF,Add=c(0,1,0,2)))
#'
#' @family makeLMEPGxData
#' @export
makeLMEPGxData <- function( x , g , chromosome = NULL , genderVar = NULL , subVar )
{  UseMethod("makeLMEPGxData",g)
}

#' @family makeLMEPGxData
#' @export
makeLMEPGxData.default<-function( x , g , chromosome = NULL , genderVar = NULL , subVar ){
  stop(paste("No method makePGxData for class:",class(g)[1])) 
}

#' @family makeLMEPGxData
#' @export
makeLMEPGxData.numeric<-function( x , g , chromosome = NULL , genderVar = NULL , subVar ){
  ids <- names( g )
  if ( min( g ) > 0 && max( g ) > 2 ){
    g <- g - min(g) # scale to 0
  }
  return( join( x[,] , eval( parse( text = paste( "data.frame(" , subVar , "= ids , Add = g )" ) ) ) , by = subVar , "left" ) )
}
#' @family makeLMEPGxData
#' @export
makeLMEPGxData.integer<-makeLMEPGxData.numeric

#' @family makeLMEPGxData
#' @export
makeLMEPGxData.character<-function( x , g , chromosome = NULL , genderVar = NULL , subVar ){
  ids <- names( g )
  Add <- makeAdditive( g , chromosome , ifelse(!is.null(genderVar),list(as.character( x[,genderVar] ) ),list(NULL))[[1]])
  return( join( x[,] , eval( parse( text = paste( "data.frame(" , subVar , "= ids , Add = Add )" ) ) ) , by = subVar , "left" ) )
}
#' @family makeLMEPGxData
#' @export
makeLMEPGxData.data.frame<-function( x , g , chromosome = NULL , genderVar = NULL , subVar ){
  g <- unlist( c( g ) )
  return( makeLMEPGxData( x , g , chromosome , genderVar , subVar ) )
}
  
#' Generate PGxG2x data
#'
#' Input an Phenotype DF, snpdata for paired patient cohort 1, anpdata for paired patient cohort 2, chromosome and gender var and output a single dataframe
#'
#' @param x data.frame of phenotypes
#' @param g vector of genotypes, In the format "#/#" or 0,1,2
#' @param g2 vector of genotypes, In the format "#/#" or 0,1,2
#' @param chromosome the chromosome the genotype is from
#' @param genderVar the colname of the gender variable in x. Make sure the values are either "Male" or "Female"
#'
#' @return a dataframe combining phenotype and genotype
#'
#' @examples
#' testDF<-data.frame(var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),gender=c("Male","Female","Male","Female"))
#' genotype1<-c("A/A","A/B","A/A","B/B")
#' genotype2<-c("A/B","A/B","A/A","A/A")
#'
#' pgx<-makePGxG2xData(testDF,genotype1,genotype2,1,"gender")
#'
#' identical(pgx,data.frame(testDF,Add=c(0,1,0,2),Add2=c(1,1,0,0)))
#'
#' @family makePGxData
#' @export
makePGxG2xData <- function( x , g , g2  , chromosome , genderVar) {
  if(length(g)!=length(g2)){
    stop("Genotypes for g and g2 are different lengths")
  }
  UseMethod("makePGxG2xData",g)
}

#' @family makePGxData
#' @export
makePGxG2xData.default<-function( x , g , g2  , chromosome , genderVar){
  stop(paste("No method makePGxData for class:",class(g)[1])) 
}

#' @family makePGxData
#' @export
makePGxG2xData.numeric<-function( x , g , g2  , chromosome , genderVar){
  if ( min( g ) > 0 && max( g ) > 2 ){
    g <- g - min(g) # scale to 0
  }
  
  PGx2<-makePGxData(rep(1,length(g2)),g2,chromosome,NULL)
  
  return( data.frame( x , Add = g, Add2 = PGx2$Add) )
}

#' @family makePGxData
#' @export
makePGxG2xData.integer<-makePGxG2xData.numeric

#' @family makePGxData
#' @export
makePGxG2xData.character<-function( x , g , g2 , chromosome , genderVar){
  Add <- makeAdditive( g , chromosome , ifelse(!is.null(genderVar),list(as.character( x[,genderVar] ) ),list(NULL))[[1]])
  PGx2<-makePGxData(rep(1,length(g2)),g2,chromosome,NULL)
  return( data.frame( x , Add = Add , Add2 = PGx2$Add) )
}

#' @family makePGxData
#' @export
makePGxG2xData.data.frame<-function( x , g , g2 , chromosome , genderVar){
  g <- unlist( c( g ) )
  return( makePGxG2xData( x , g , g2 , chromosome , genderVar ) )
}

#' Generate LME PGxG2x data
#'
#' Input an Phenotype DF, snpdata for paired patient cohort 1, anpdata for paired patient cohort 2, chromosome and gender var and output a single dataframe
#'
#' @param x data.frame of phenotypes
#' @param g vector of genotypes, In the format "#/#" or 0,1,2
#' @param g2 vector of genotypes, In the format "#/#" or 0,1,2
#' @param chromosome the chromosome the genotype is from
#' @param genderVar the colname of the gender variable in x. Make sure the values are either "Male" or "Female"
#'
#' @return a dataframe combining phenotype and genotype
#'
#' @examples
#' testDF<-data.frame(ID=paste0("sub",1:4),var1=c(1,1,1,1),var2=c(2,2,2,2),var3=c(3,3,3,3),gender=c("Male","Female","Male","Female"))
#' genotype1<-c(sub1="A/A",sub2="A/B",sub3="A/A",sub4="B/B")
#' genotype2<-c("A/B","A/B","A/A","A/A")
#' 
#' pgx<-makeLMEPGxPGx2Data(testDF,genotype1,genotype2,1,"gender","ID")
#'
#' identical(pgx,data.frame(testDF,Add=c(0,1,0,2),Add2=c(1,1,0,0)))
#'
#' @family makeLMEPGxData
#' @export
makeLMEPGxPGx2Data <- function( x , g , g2 , chromosome = NULL , genderVar = NULL , subVar ){
  if(length(g)!=length(g2)){
    stop(paste0("Genotypes for g and g2 are different lengths: ",length(g),", ",length(g2)))
  }
  UseMethod("makeLMEPGxPGx2Data",g)
}

#' @family makeLMEPGxData
#' @export
makeLMEPGxPGx2Data.default<-function( x , g , g2 , chromosome = NULL , genderVar = NULL , subVar ){
  stop(paste("No method makePGxData for class:",class(g)[1])) 
}

#' @family makeLMEPGxData
#' @export
makeLMEPGxPGx2Data.numeric<-function( x , g , g2 , chromosome = NULL , genderVar = NULL , subVar ){
  
  ids <- names( g )
  if ( min( g ) > 0 && max( g ) > 2 ){
    g <- g - min(g) # scale to 0
  }
  
  PGx2<-makePGxData(rep(1,length(g2)),g2,chromosome,NULL)
  
  return( join( x[,] , eval( parse( text = paste( "data.frame(" , subVar , "= ids , Add = g , Add2 = PGx2$Add )" ) ) ) , by = subVar , "left" ) )
}
#' @family makeLMEPGxData
#' @export
makeLMEPGxPGx2Data.integer<-makeLMEPGxPGx2Data.numeric

#' @family makeLMEPGxData
#' @export
makeLMEPGxPGx2Data.character<-function( x , g , g2 , chromosome = NULL , genderVar = NULL , subVar ){
  ids <- names( g )
  Add <- makeAdditive( g , chromosome , ifelse(!is.null(genderVar),list(as.character( x[,genderVar] ) ),list(NULL))[[1]])
  
  PGx2<-makePGxData(rep(1,length(g2)),g2,chromosome,NULL)
  
  return( join( x[,] , eval( parse( text = paste( "data.frame(" , subVar , "= ids , Add = Add , Add2 = PGx2$Add )" ) ) ) , by = subVar , "left" ) )
}

#' @family makeLMEPGxData
#' @export
makeLMEPGxPGx2Data.data.frame<-function( x , g , g2 , chromosome = NULL , genderVar = NULL , subVar ){
  g <- unlist( c( g ) )
  return( makeLMEPGxPGx2Data( x , g , g2 , chromosome , genderVar , subVar ) )
}

#' Generate genotype calls for PGx data
#'
#' Input an Genotype calls, Chromosome and gender, returns genotype calls in numeric format
#'
#' @param g vector of genotypes, In the format "#/#",
#' @param chromosome the chromosome the genotype is from
#' @param gender the vector of gender. Make sure the values are either "Male" or "Female"
#'
#'
#' @return a vector of genotype calls numerically
#'
#' @examples
#' gender<-c("Male","Female","Male","Female")
#' genotype<-c("A/A","A/B","A/A","B/B")
#'
#' add<-makeAdditive(genotype,1,gender)
#'
#' identical(add,c(0,1,0,2))
#'
#' @export
makeAdditive <- function( g , chromosome , gender )
{
  alleles <- table( unlist( lapply( as.character( g ) , function( x ) strsplit( x , "/" , fixed = TRUE ) ) ) )
  if(names(which.max(alleles))!=names(which.min(alleles))){
  genotypeClasses <- c( paste( names( alleles )[which.max( alleles )] , names( alleles )[which.max( alleles )] , sep = "/" ) ,
                        paste( names( alleles )[which.max( alleles )] , names( alleles )[which.min( alleles )] , sep = "/" ) ,
                        paste( names( alleles )[which.min( alleles )] , names( alleles )[which.max( alleles )] , sep = "/" ) ,
                        paste( names( alleles )[which.min( alleles )] , names( alleles )[which.min( alleles )] , sep = "/" ) )
  genotypeClasses<-genotypeClasses[c(1,which(genotypeClasses[c(2,3)]%in%levels(factor(as.character(g))))+1,4)]
  }else{
    genotypeClasses<-levels(factor(as.character(g)))
  }
  
  if(is.null(gender))
    gender<-"Female"
  
  return( unlist( apply( cbind( g , gender ) , 1 , additiveGenotype , genotypeClasses , chromosome ) ) )
}
#' Generate genotype calls for makeAdditive function
#'
#' Input an Genotype call, Chromosome and gender, returns a genotype call in numeric format
#'
#' @param x vector of two elements, [1]genotype call and [2]gender
#' @param genotypeClasses the different levels that the genotypes can manifest in
#' @param chromosome the chromosome of the genotype
#'
#'
#' @return a list, with a model.frame object and a data.frame
#'
#' @examples
#' gender<-c("Male","Female","Male","Female")
#' genotype<-c("A/A","A/B","A/A","B/B")
#'
#' genotypeClasses<-c("A/A","A/B","B/B")
#'
#' add<-additiveGenotype(cbind(genotype,gender)[1,],genotypeClasses,1)
#' identical(add,0)
#' 
#' add2<-apply( cbind( genotype , gender ) , 1 , additiveGenotype , genotypeClasses , 1 ) 
#' identical(add2,c(0,1,0,2))
#'
#' @export
additiveGenotype <- function( x , genotypeClasses , chromosome )
{
  if ( is.na( x[1] ) )
  {
    return( NA )
  }
  else
  {
    if ( chromosome == "X" )
    {
      geneFrame <- data.frame( Genotype = genotypeClasses , Female = c( 0 , 1 , 2 ) , Male = c( 0 , NA , 2 ) )
      return( ifelse( x[2] == "Female" , geneFrame[which( x[1] == geneFrame[,1] ),2] , geneFrame[which( x[1] == geneFrame[,1] ),3] ) )
    }
    else
    {
      geneFrame <- data.frame( Genotype = genotypeClasses , Additive = c( 0 , 1 , 2 ) )
      return( geneFrame[which( x[1] == geneFrame[,1] ),2] )
    }
  }
}


#' Run a GWAS
#'
#' Run the GWAS using the function identified in func. Can be ran in parallel or series
#'
#' @param pgx DF of phenotypes for each subject
#' @param snp Matrix of genotype calls, snp by row, subject by column. Genotype format varies by called "func".
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param func The name of the function to apply in the GWAS. This is meant to be a user-defined function, with the mandatory output of a list containing at least one value called "pval". Inluded in the Package are some examples.
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param adjust The method to adjust the P-values by after the GWAS. defaults to "BY", see \code{\link[stats]{p.adjust}} for more options
#' @param replace The name of the variable to remove from the Model if necessary
#' @param snowFall Boolean Value, Should the GWAS be run in parallel or series. If TRUE, use \code{\link[snowfall]{sfInit}} prior to \code{runGWAS}
#' @param file Name of file which to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param subCoeff coefficients to keep and to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param subEffect Effect to keep and to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param dbPath path to database to save results in if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param dbConnection a sqlite connection to a DB, alternative to dbPath to save results if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param dbTABLENAME The name of the table in the DB to save the results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param ... additional arguments passed to func 
#'
#' @return a list, with a data.frame of P-Values, and a list of models for each snp
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasRes<-runGWAS(pgx=examplePGx,snp=exampleSNP,Model=exampleModel,func = gwasFunc,interactionTerm="var2",genderVar="gender",databasePath="runGWAS.DB", dbTable="runGWAS_Test", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
#' @export
runGWAS <- function( pgx , snp , Model , func = gwasFunc , interactionTerm = "ARMACD" , genderVar = "SEX.Decode" , adjust = "BY" , replace = NULL , snowFall = FALSE , file = NULL , subCoef = NULL , subEffect = NULL , dbPath = NULL , dbConnection = NULL , dbTABLENAME = NULL , overwriteEntry = NULL , tempFile = tempfile() , ... )# family = NULL , randomForm = NULL , correlationForm = NULL , subjectID = NULL , ... )
{
  if ( length( interactionTerm ) < 1 )
  {
    interactionTerm <- NULL
  }
  if ( !is.null( replace ) )
  {
    Model <- model.frame( updateModel( pgx , Model , NULL , replace ) , pgx )
  }
  if ( inherits( snp , c( "CNSet" , "SnpSuperSet" ) ) )
  {
    fd <- calls( snp )
    finalizer( fd ) <- "close"
    if ( snowFall )
    {
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "ff" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( ff , keep.source = FALSE )
      }
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "car" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( car , keep.source = FALSE )
      }
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "nnet" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( nnet , keep.source = FALSE )
      }
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "lme4" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( lme4 , keep.source = FALSE )
      }
      #Add AxioSerializer and RSQLite
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "AxioSerializer" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( AxioSerializer , keep.source = FALSE )
      }
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "RSQLite" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( RSQLite , keep.source = FALSE )
      }
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "survival" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( survival , keep.source = FALSE )
      }
      
      sfExport( "fd" )
      sfClusterEval( open( fd ) )
      sfExport( as.character( quote( gwasFunc ) ) )
      sfExport( paste( as.character( quote( gwasFunc ) ) , "sf" , sep = "" ) )
      sfExport( "gwasFuncsf" )
      sfExport( "gwasLMEFuncsf" )
      sfExport( "gwasLMEHFuncsf" )
      sfExport( "gwasGLMFuncsf" )
      sfExport( "gwasMGLMFuncsf" )
      sfExport( "gwasCOXPHFuncsf" )
      sfExport( "makePGxData" )
      sfExport( "makeAdditive" )
      sfExport( "additiveGenotype" )
      sfExport( "updateModel" )
      sfExport( "reducedModel" )
      sfExport( "writeResultsSNPWise" )
      sfExport( "getCoefficients" )
      sfExport( "makeCoefLabels" )
      sfExport( "alleleSplit" )
      sfExport( "resultsDataFrame" )
      sfExport( "mafFreqIllumina" )
      sfExport( "hardyTest" )
      sfExport( "hardyTestIllumina" )
      sfExport( "hardyTestIlluminaSF" )
      tpgx <- data.frame( pgx , Add = runif( nrow( pgx ) , 0 , 2 ) )
      depVar <- rownames( attr( terms( Model ) , "factors" ) )[attr( terms( Model ) , "response" )]
      if ( inherits( tpgx[,depVar] , "factor" ) )
      {
        if ( length( levels( tpgx[,depVar] ) ) > 2 )
        {
          levs <- levels( tpgx[,depVar] )[-1]
        }
        else
        {
          levs <- NULL
        }
      }
      else
      {
        levs <- NULL
      }
      cfnts <- makeCoefLabels( attr( terms( formula( updateModel( tpgx , Model , interactionTerm ) ) ) , "term.labels" ) , tpgx , levs )
      if ( !is.null( subEffect ) )
      {
        indx <- c()
        for ( i in subEffect )
        {
          indx <- c( indx , grep( i , cfnts , fixed = TRUE ) )
        }
        cfnts <- cfnts[indx]
      }
      dfTmp <- as.data.frame( matrix( nrow = 1 , ncol = ( 3 * length( cfnts ) + 4 ) ) )
      names( dfTmp ) <- c( "HWE" , "MAF" , "Adjusted P-Value" , "Raw P-Value" , paste( "Estimate" , cfnts ) , paste( "Std. Err." , cfnts ) , paste( "p-value" , cfnts ) )
      if ( !is.null( subCoef ) )
      {
        indx <- c( seq( 1 , 4 ) )
        for ( i in subCoef )
        {
          indx <- c( indx , grep( i , names( dfTmp ) , fixed = TRUE ) )
        }
        dfTmp <- dfTmp[,indx]
      }
#       ffCall <- paste( "ffObject <- ffdf( row.names = rownames( fd ) , " , paste( dQuote( names( dfTmp ) ) , collapse = " = ff( vmode = \"double\" , length = nrow( fd ) ) , " ) , "= ff( vmode = \"double\" , length = nrow( fd ) ) )" )
#       eval( parse( text = ffCall ) )
#       sfExport( "ffObject" )
#       sfClusterEval( open( ffObject ) )
#       gwasRes <- sfLapply( chunk.ffdf( fd , by = nrow( fd ) / length( sfGetCluster() ) ) , func , fd , pgx , Model , interactionTerm , subCoef , subEffect , ffObject , ... )
#       sfClusterEval( close( fd ) )
#       ffObject[,"Adjusted.P.Value"] <- p.adjust( ffObject[,"Raw.P.Value"] , adjust )
# #      ffObject[,"HWE"] <- p.adjust( ffObject[,"HWE"] , adjust )
# #      ffObject <- ffObject[ffdforder( ffObject[4] ),]
# #      print( class( ffObject ) )
#       write.csv( file = file , ffObject )
#       sfClusterEval( close( ffObject ) )
      
      #############
      #### Create SQLite DB connection
      #############
      if( !is.null(dbConnection) && inherits(dbConnection,c("SQLiteConnection","RSQLite"))){
        if(is.null(dbPath)){
          error("Please Provide either the path to a database or the RSQLite db Connection")
        }else{
          dbConnection<-dbConnect(SQLite(),dbPath)
        }
      }
      
      #############
      #### Create Table to store stats of Analysis
      #############
      if(is.null(dbTABLENAME)){
        dbTABLENAME<-"SNP_STATS"
      }
      if(dbTABLENAME%in%dbListTables(dbConnection)){
        dbSendQuery(dbConnection,paste("DROP TABLE",dbTABLENAME))
        dbSendQuery(dbConnection,paste("DROP TABLE",paste( "Anova" , dbTABLENAME , sep = "_" )))
      }
      
      #Create a table with all the columns in dfTmp by Adjusted P-Value to allow merging later
      databaseQuery<-paste("CREATE TABLE",dbTABLENAME,"(SNP varchar ,",paste0('"',setdiff(names(dfTmp),"Adjusted P-Value"),'"',collapse=" double , "),");")
      dbSendQuery(dbConnection,databaseQuery)
      databaseQuery<-paste("CREATE TABLE",paste( "Anova" , dbTABLENAME , sep = "_" ),"(SNP varchar ,",paste0('"',setdiff(names(dfTmp),"Adjusted P-Value"),'"',collapse=" double , "),");")
      dbSendQuery(dbConnection,databaseQuery)
      dbDisconnect(dbConnection)
      
      #export function to be sfLappy-ed to each worker
      eval(parse(text=paste("sfExport('",func,"')")))
      
      gwasRes <- sfLapply( chunk.ffdf( fd , by = nrow( fd ) / length( sfGetCluster() ) ) , func , fd , pgx , Model , interactionTerm , subCoef , subEffect , dbPath , dbTABLENAME , ... )
      
      #reconnect to write 
      dbConnection<-dbConnect(SQLite(),dbConnection@dbname)
      
      #get Adjusted pvalue
      pvalDB<-dbGetQuery(dbConnection,paste0("SELECT SNP,[Raw P-Value] FROM ",dbTABLENAME))
      pvalDB$`Adjusted P-Value` <-  p.adjust( pvalDB$`Raw P-Value`, adjust)
      
      dbWriteTable(dbConnection,paste0(dbTABLENAME,"_ADJPVAL"),pvalDB)
      
      #Join Tables into a merged table, and then filter out duplicated Columns, remove extra tables to become final table
      dbGetQuery(dbConnection,paste0("CREATE TABLE ",paste0(dbTABLENAME,"_MERGED")," AS SELECT * FROM ",dbTABLENAME," INNER JOIN ",paste0(dbTABLENAME,"_ADJPVAL")," ON ",dbTABLENAME,".SNP=",paste0(dbTABLENAME,"_ADJPVAL"),".SNP"))
      names<-names(dbGetQuery(dbConnection,paste("SELECT * FROM",paste0(dbTABLENAME,"_MERGED"),"WHERE 1==0")))
      newTableQuery<-paste0("CREATE TABLE ",paste0(dbTABLENAME,"_FINAL")," AS SELECT [",paste(names[-grep(":1",names)],collapse="],["),"] FROM ",paste0(dbTABLENAME,"_MERGED"))
      dbGetQuery(dbConnection,newTableQuery)
      lapply(c(paste0(dbTABLENAME,c("","_ADJPVAL","_MERGED"))),function(x,dbcon)dbGetQuery(dbcon,paste("DROP TABLE",x)),dbConnection)
      dbGetQuery(dbConnection,paste("ALTER TABLE",paste0(dbTABLENAME,"_FINAL"),"RENAME TO",dbTABLENAME))
      dbGetQuery(dbConnection,paste("DROP TABLE",paste0(dbTABLENAME,"_FINAL")))
      
      
      modelFits<-dbReadTable(dbConnection,dbTABLENAME)
      rownames(modelFits)<-modelFits$SNP
      adjustedPValue<-dbGetQuery(dbConnection,paste("SELECT SNP,[Adjusted P-Value] FROM ",dbTABLENAME))
      rownames(adjustedPValue)<-adjustedPValue$SNP
      adjustedPValue<-as.numeric(adjustedPValue$`Adjusted P-Value`,names=adjustedPValue$SNP)
      
      dbDisconnect(dbConnection)
      
      write.csv( file = file , ffObject )
      
      return( list( adjPVals = ffObject[,"Adjusted.P.Value"] , modelFits = ffObject ) )
    }
    else
    {
      gwasRes <- lapply( seq( 1 , dim( snp )["Features"] ) , func , fd , pgx , Model , interactionTerm , genderVar , ... )
      close( fd )
      adjPVals <-  p.adjust( unlist( lapply( gwasRes , function( x ) x$pVal ) ) , adjust )
      return( list( adjPVals = adjPVals , modelFits = gwasRes ) )
    }
  }
  else
  {
    if ( snowFall )
    {
      #Add AxioSerializer and RSQLite
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "AxioSerializer" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( AxioSerializer , keep.source = FALSE )
      }
      lidx <- split( seq( 1 , nrow( snp ) ) , cut( seq_along( seq( 1 , nrow( snp ) ) ) , sfNodes() , labels = FALSE ) )
      gwasFunc <- "gwasRes <- sfLapply( lidx , func , snp , pgx , Model , interactionTerm , genderVar , dbPath , dbTABLENAME , overwriteEntry , tempFile"

    }
    else
    {
      gwasFunc <- "gwasRes <- lapply( seq( 1 , nrow( snp ) ) , func , snp , pgx , Model , interactionTerm , genderVar, dbPath , dbTABLENAME , overwriteEntry , tempFile"

    }
    
    #names( gwasRes ) <- rownames( snp )
    
    Args <- names( formals( func ) )[-c( 1:10 )]
    allArgs <- names( as.list( sys.call() ) )
    missingArgs <- setdiff( Args , allArgs )
    functionArgs <- list( ... )[Args]

    if ( length( missingArgs ) > 0 )
    {
      for ( i in seq_along( missingArgs ) )
      {
        print( paste( "Argument" , missingArgs[i] , "is missing with no default in apply function." ) )
      }
      stop( "Missing arguments in function." )
    }
    if ( length( Args ) > 0 )
    {
      for ( i in seq_along( Args ) )
      {
        gwasFunc <- paste( gwasFunc , "," , paste( Args[i] , ifelse( inherits( functionArgs[[Args[i]]] , "corStruct" ) , getCor( functionArgs[[Args[i]]] ) ,
                                                                     ifelse( inherits( functionArgs[[Args[i]]] , c( "formula" , "call" ) ) , deparse( functionArgs[[Args[[i]]]] ) ,
                                                                             ifelse( inherits( functionArgs[[Args[i]]] , c("data.frame","ff","dbDF","matrix","data.table" ) ) , paste0("functionArgs[['",Args[[i]],"']]") ,
                                                                                     paste( "\"" , functionArgs[[Args[i]]] , "\"" , sep = "" ) ) ) ) , sep = " = "  ) )

      }
    }
    gwasFunc <- paste( gwasFunc , ")" )
    
    eval( parse( text = gwasFunc ) )
    
    adjPVals <-
      data.frame(do.call('rbind' ,
        lapply(gwasRes , function(x1){
            if (is.list(x1)) {  do.call('rbind', lapply(x1, function(x2)
                if (is.list(x2)) { do.call('c', x2) } else{ x2 }))
            } else{ x1 }}
            
        )))
        
    
    adjPVals$pVal<-as.numeric(as.character(adjPVals$pVal))
    adjPVals$adjpval <- p.adjust( as.numeric( as.character( adjPVals$pVal ) ) , adjust )

    rownames(adjPVals)<-adjPVals$ID
    colnames(adjPVals)<-c("pVal","SNPID","adjpval")
    

    # return( list( adjPVals = adjPVals , modelFits = gwasRes ) )
    return( adjPVals )
    
  }
}

#' Clean GWAS results and store in table
#'
#' A helper function called by gwas_XXX_FuncSF in runGWAS to extract data and store results in a db
#'
#' @param x list of lists generated from gwas_XXX_Funcsf helper functions
#' @param dbPath The path to database to save results
#' @param TABLE The table in the db to save the results to
#' 
#' @return NA
#'
#' @import stats
#' @import AxioSerializer
#' @import car
returnDF_writeDB<-function(x, dbPath, TABLE){
  
  dbConnection <- RSQLite::dbConnect(RSQLite::SQLite(), dbPath)
  
  dataFrame <- as.data.frame(matrix(unlist(lapply(x , function(x, dbConn, dbTable) {
    #Code to insert into table
    databaseQuery <- paste("INSERT INTO",dbTable,"(SNP,",paste0('"', c(names(x$target),names(x$data)), '"', collapse = ","),") VALUES(",paste(c(x$target, x$data), collapse = ","),");")
    databaseQuery <- gsub(",NA,", ",NULL,", databaseQuery)
    repeat {
      rv <- try(RSQLite::dbGetQuery(dbConn, databaseQuery))
      if (!is(rv, "try-error"))
        break
    }
    return(as.matrix(x$data))
  },
  dbConnection,
  TABLE)),
  nrow = length(pVals),
  byrow = TRUE),
  stringsAsFactors = FALSE)
  
  rownames(dataFrame) <- names(x$pVal)  
  
  return()
}
#' Clean LM objects of all kinds results and store in table
#'
#' A helper function to reduce the size of the lm objects. They had connections to the environment the lm object was created in, and would become large when serialized
#'
#' @param LMINPUT any sort of lm object
#' 
#' @return stripped down lm object
#' @export
stripGLM_callEnv <- function(LMINPUT){
  # for some reason, the LM object was storing the call environment, which could get rather LARGE, 
  # if the input pgx was large, which since it has to do with SNPS, often was
  attributes(LMINPUT$terms)$'.Environment'<-c()
  attributes(attributes(LMINPUT$model)$terms)$'.Environment'<-c()
  LMINPUT$data = c()
  return(LMINPUT)
}

#' Clean LME objects of all kinds results and store in table
#'
#' A helper function to reduce the size of the lme objects. They had connections to the environment the lme object was created in, and would become large when serialized
#'
#' @param LMEINPUT lme object
#' 
#' @return stripped down lm object
#' @export
stripLME_callEnv <- function(LMEINPUT){
  # for some reason, the LM object was storing the call environment, which could get rather LARGE, 
  # if the input pgx was large, which since it has to do with SNPS, often was
  attributes(LMEINPUT$terms)$'.Environment'<-c()
  attributes(attributes(LMEINPUT$model)$terms)$'.Environment'<-c()
  attributes(attributes(LMEINPUT$model$corStruct)$formula)$.Environment<-NULL
  attributes(attributes(LMEINPUT$model$reStruct$STUDYID)$formula)$.Environment<-NULL
  attributes(attributes(LMEINPUT$modelStruct$corStruct)$formula)$.Environment<-NULL
  attributes(attributes(LMEINPUT$modelStruct$reStruct$STUDYID)$formula)$.Environment<-NULL
  LMEINPUT$data = c()
  return(LMEINPUT)
}

#' Clean COXPH objects of all kinds results and store in table
#'
#' A helper function to reduce the size of the coxph objects. They had connections to the environment the coxph object was created in, and would become large when serialized
#'
#' @param COXPHINPUT coxph object
#' 
#' @return stripped down lm object
#' @export
stripCOXPH_callEnv <- function(COXPHINPUT){
  # for some reason, the LM object was storing the call environment, which could get rather LARGE, 
  # if the input pgx was large, which since it has to do with SNPS, often was
  attributes(COXPHINPUT$terms)$'.Environment'<-c()
  attributes(COXPHINPUT$formula)$'.Environment'<-c()
  COXPHINPUT$data = c()
  return(COXPHINPUT)
}

#' Clean LME objects of all kinds results and store in table
#'
#' A helper function to reduce the size of the lme objects. They had connections to the environment the lme object was created in, and would become large when serialized
#'
#' @param LMEINPUT lme object
#' 
#' @return stripped down lm object
#' @export
stripGLME_callEnv <- function(GLMEINPUT){
  # for some reason, the LM object was storing the call environment, which could get rather LARGE, 
  # if the input pgx was large, which since it has to do with SNPS, often was
  GLMEINPUT@resp@.xData<-new.env()
  attributes(attributes(GLMEINPUT@frame)$formula)$.Environment<-NULL

  return(GLMEINPUT)
}