#' run a Linear Model for a GWAS
#'
#' A function called by runGWAS to perform linear models for each snp
#'
#' @param i what row from SNP to call the genotypes from
#' @param g Matrix of "#/#" genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param dbPath The path to database to save results
#' @param dbTable The name of the table in the DB to save the results to
#' @param snpInfo Details about each snp. Requires at least a column titled 'Chr', for chromosome of snp
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a list with the P-Value, and a list of models for each snp
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasFuncRes<-gwasFunc(i=1,g=exampleSNP, pgx=examplePGx,Model=exampleModel,interactionTerm="var2",genderVar="gender",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
#'
#' @export
gwasFunc <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , snpInfo )
{
  pgx <- makePGxData( pgx , g[i,] , snpInfo[i,"Chr"] , genderVar )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- lm( newForm , data = pgx )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "LM_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
    writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "LM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>F)"[which( rownames( aTable ) == "Add" )], ID=rownames( g )[i]  ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>F)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )], ID=rownames( g )[i]  ) )
  }
}

#' run a LM for a GWAS in Parallel
#'
#' A function called by runGWAS to perform LM for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param subCoef The name of the sub Coeffects varable to include in results
#' @param subEffect The name of the sub Effects varable to include in results
#' @param dbPath The path to database to save results
#' @param Table The name of the table in the DB to save the results to
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a vector of P values, and saves lm models to the table/db specified in DBPath/Table
#'
#' @import stats
#' @export
gwasFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  return( lapply( i , gwasFuncsf , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = overwriteEntry , tempFile = tempFile ) )
}

#' run a LM for a GWAS
#'
#' A helper function called by gwasFuncSF in runGWAS to perform a LM for each snp
#'
#' @param i index of snp of interest in g
#' @param g matrix object of numeric genotype calls, snp by row, sample by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param subCoef The name of the sub Coeffects varable to include in results
#' @param subEffect The name of the sub Effects varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 3 entries, a P value, dataframe of results and the name of the target. LM model is also saved to the db specified in dbPath
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasFuncsf <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- makePGxData( pgx , g[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- lm( newForm , data = pgx )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "LM_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
    writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "LM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>F)"[which( rownames( aTable ) == "Add" )], ID=rownames( g )[i]  ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>F)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )], ID=rownames( g )[i] ) )
  }
}

#' run a lm for a GWAS in Parallel or series
#'
#' A function called by runGWAS to perform a lm for each snp, designed for both SnowFall parallel processing and series.
#'
#' @param i the index of the snp in g to run the lm on
#' @param g matrix object of numeric genotype calls, snp by row, sample by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the genderVar in pgx
#' @param dbPath The path to database to save results
#' @param dbTable the name of the table to save the lm objects to. done using \code{\link[AxioSerializer]{writeObjectToTable}}
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with two entries, a a P value, and lm model. lm model is also saved to the db specified in databasePath.
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasFuncFF <- function( i , g , pgx , Model , interactionTerm , genderVar, dbPath=NULL, dbTable=NULL , overwriteEntry=FALSE , tempFile = tempfile() )
{
  Add<-g[i,]-1
  if(is.data.frame(Add))
    Add<-do.call('c',Add)
  pgx <- data.frame( pgx , Add = Add )
  
  newForm <- updateModel( pgx , Model , interactionTerm , genderVar )
  lmTmp <- lm( newForm , data = pgx )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "LM_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
    writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "LM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>F)"[which( rownames( aTable ) == "Add" )], ID=rownames( g )[i]  ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>F)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )], ID=rownames( g )[i] ) )
  }
}

#' run a Linear Mixed Effects Model for a GWAS
#'
#' A function called by runGWAS to perform Linear Mixed Effects models for each snp
#'
#' @param i what row from SNP to call the genotypes from
#' @param g Matrix of "#/#" genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param randomTerm The name of the random Term in the pgx data.frame
#' @param heterpTerm The name of the hetero Term in the pgx data.frame
#' @param dbPath The path to database to save results
#' @param dbTable The name of the table in the DB to save the results to
#' @param snpInfo Details about each snp. Requires at least a column titled 'Chr', for chromosome of snp
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P-Value or an NA if unable to calculate
#'
#' @import nlme
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),var4=c(1,2,1,1,2,1,1,2),var5=c(2,2,2,1,1,1,2,1),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasLMEHFuncRes<-gwasLMEHFunc(i=1,g=exampleSNP, pgx=examplePGx,Model=exampleModel,interactionTerm="var2",genderVar="gender",randomTerm="var4",heteroTerm="var5",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
gwasLMEHFunc <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , randomTerm , heteroTerm , snpInfo )
{
  pgx <- makePGxData( pgx , g[i,] , snpInfo[i,"Chr"] , genderVar )
  newForm <- formula( updateModel( pgx , Model , interactionTerm ) )
  randForm <- eval( parse( text = paste( "list(" , randomTerm , "= ~ 1 |" , heteroTerm , ")" ) ) )
  hetForm <- eval( parse( text = paste( "~ 1 |" , heteroTerm ) ) )
  lmTmp <- try( nlme::lme( newForm , data = pgx , random = randForm , weights = nlme::varIdent( form = hetForm ) , na.action = "na.exclude" ) )
  lmWrite <- ifelse( inherits( lmTmp , "try-error" ) , NA , stripGLM_callEnv( lmTmp ) )
  aTable <- try( car::Anova(  lmTmp , type = "III" , singular.ok = TRUE ) )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( lmTmp , paste0( "LMEH_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
    writeObjectToTable( aTable , paste0( "LMEH_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  if ( inherits( lmTmp , "try-error" ) )
  {
    return( list( pVal = NA ) )
  }
  else
  {
    if ( is.null( interactionTerm ) )
    {
      return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( lmTmp ) == "Add" )], ID=rownames( g )[i]  ) )
    }
    else
    {
      return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( lmTmp ) == paste( interactionTerm , "Add" , sep = ":" ) )], ID=rownames( g )[i]  ) )
    }
  }
}

#' run a LMEH for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LMEH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomTerm The name of the random varable to include in results
#' @param heterogeneousTerm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a vector of P values, and saves lm models to the db specified in DBPath
#' @export
gwasLMEHFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() , randomTerm , heterogeneousTerm )
{
  return( sapply( i , gwasLMEHFuncsf , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry , tempFile , randomTerm , heterogeneousTerm ) )
}

#' run a LMEH for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LMEH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomTerm The name of the random varable to include in results
#' @param heterogeneousTerm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P value, and saves lme models to the db specified in dbPath
#'
#' @import nlme
#' @import AxioSerializer
gwasLMEHFuncsf <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() , randomTerm , heterogeneousTerm )
{
  pgx <- makePGxData( pgx , g[i,] , snpInfo[i,"Chr"] , genderVar )
  newForm <- formula( updateModel( pgx , Model , interactionTerm ) )
  randForm <- eval( parse( text = paste( "list(" , randomTerm , "= ~ 1 |" , heteroTerm , ")" ) ) )
  hetForm <- eval( parse( text = paste( "~ 1 |" , heteroTerm ) ) )
  lmTmp <- try( nlme::lme( newForm , data = pgx , random = randForm , weights = nlme::varIdent( form = hetForm ) , na.action = "na.exclude" ) )
  lmWrite <- ifelse( inherits( lmTmp , "try-error" ) , NA , stripGLM_callEnv( lmTmp ) )
  aTable <- try( car::Anova(  lmTmp , type = "III" , singular.ok = TRUE ) )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( lmTmp , paste0( "LMEH_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
    writeObjectToTable( aTable , paste0( "LMEH_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  if ( inherits( lmTmp , "try-error" ) )
  {
    return( list( pVal = NA ) )
  }
  else
  {
    if ( is.null( interactionTerm ) )
    {
      return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( lmTmp ) == "Add" )], ID=rownames( g )[i]  ) )
    }
    else
    {
      return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( lmTmp ) == paste( interactionTerm , "Add" , sep = ":" ) )], ID=rownames( g )[i]  ) )
    }
  }
}

#' run a Linear Mixed Effects Model for a GWAS
#'
#' A function called by runGWAS to perform Linear Mixed Effects models for each snp
#'
#' @param i what row from SNP to call the genotypes from
#' @param g Matrix of "#/#" genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param randomTerm The name of the random Term in the pgx data.frame
#' @param heterpTerm The name of the hetero Term in the pgx data.frame
#' @param dbPath The path to database to save results
#' @param dbTable The name of the table in the DB to save the results to
#' @param snpInfo Details about each snp. Requires at least a column titled 'Chr', for chromosome of snp
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P-Value or an NA if unable to calculate
#'
#' @import nlme

#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),var4=c(1,2,1,1,2,1,1,2),var5=c(2,2,2,1,1,1,2,1),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasLMEHFuncRes<-gwasLMEHFunc(i=1,g=exampleSNP, pgx=examplePGx,Model=exampleModel,interactionTerm="var2",genderVar="gender",randomTerm="var4",correlationTerm="var5",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
gwasLMEFunc <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , randomForm , corrFun , phiVal , corrForm , snpInfo )
{
  pgx <- AxioGWAS:::makeLMEPGxData( pgx , g[i,] , NULL , NULL , subVar = subjectID )
  correlationForm <- parse( text = AxioGWAS:::getCorList( list( corrForm = corrForm , phiVal = phiVal , corrFun = corrFun ) ) )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- eval( parse( text = paste( "lme( fixed =" , newForm , ", data = pgx , random =" , randomForm , ", correlation =" , correlationForm , ", na.action = na.exclude )" ) ) )
#  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  aTable <- anova( lmTmp )
  writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "LME_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( stripGLM_callEnv( car::Anova( lmTmp , type = "III" , singular.ok = TRUE ) ) , paste0( "LME_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
#    return( list(  pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , data = aTable ) )
    return( list(  pVal = aTable$"p-value"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i]) )
  }
  else
  {
#    return( list(  pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , data = aTable ) )
    return( list(  pVal = aTable$"p-value"[which( rownames( aTable ) == "Add" )], ID=rownames( g )[i]) ) 
  }
}

#' run a LME for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LME for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomForm The name of the random varable to include in results
#' @param correlationForm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a vector of P values, and saves lm models to the db specified in DBPath
#' @export
gwasLMEFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , randomForm , corrFun , phiVal , corrForm , subjectID )
{
  return( lapply( i , gwasLMEFuncsf , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry , tempFile , randomForm , corrFun , phiVal , corrForm , subjectID ) )
}

#' run a LME for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LMEH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomForm The name of the random varable to include in results
#' @param correlationForm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P value, and saves lme models to the db specified in dbPath
#'
#' @import nlme
#' @import AxioSerializer
gwasLMEFuncsf <- function( i , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() , randomForm , corrFun , phiVal , corrForm , subjectID ){
  pgx <- AxioGWAS:::makeLMEPGxData( pgx , g[i,] , NULL , NULL , subVar = subjectID )
  correlationForm <- parse( text = AxioGWAS:::getCorList( list( corrForm = corrForm , phiVal = phiVal , corrFun = corrFun ) ) )
  newForm <- updateModel( pgx , Model , interactionTerm )
  #wrapped with a "try"
  lmTmp <- try(eval( parse( text = paste( "lme( fixed = newForm , data = pgx , random =" , randomForm , ", correlation =" , correlationForm , ", na.action = na.exclude )" ) ) ) )
  
  # if it fails, write why it failed for both regression and anova objects, set pval returned to "NA"
  if(inherits(lmTmp,"try-error")){
    lmTmp<-as.character(lmTmp)
    aTable<-lmTmp
    pval<-NA
  }else{
    aTable <- anova( lmTmp )
    lmTmp<-stripLME_callEnv( lmTmp )
    if ( is.null( interactionTerm ) ){
      pval<-aTable$"p-value"[which( rownames( aTable ) == "Add" )]
    }else{
      pval<-aTable$"p-value"[which( rownames( aTable ) ==  paste( interactionTerm , "Add" , sep = ":" ) )]
    }
  }
  
  #save values to tables
  writeObjectToTable( lmTmp , paste0( "LME_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( aTable , paste0( "LME_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  
  return( list(  pVal = pval , ID=rownames( g )[i] ) )
}


#' run a LME for a GWAS in Parallel with two input genotypes
#'
#' A function called by runGWAS to perform a LME for each snp, comparing two genotypes, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomForm The name of the random varable to include in results
#' @param correlationForm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#' @param g2 object of genotype calls, snp by row, subject by column, for the second genotype to compare
#' 
#' @return a vector of P values, and saves lm models to the db specified in DBPath
#' @export
gwasLMEFuncSF_2Genotypes <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , randomForm , corrFun , phiVal , corrForm , subjectID , g2)
{
  return( lapply( i , gwasLMEFuncsf_2Genotypes , g , g2 , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry , tempFile , randomForm , corrFun , phiVal , corrForm , subjectID ) )
}

#' run a LME for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LMEH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomForm The name of the random varable to include in results
#' @param correlationForm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#' @param g2 object of genotype calls, snp by row, subject by column, for the second genotype to compare

#' @return a P value, and saves lme models to the db specified in dbPath
#'
#' @import nlme
#' @import AxioSerializer
gwasLMEFuncsf_2Genotypes <- function( i , g , g2, pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() , randomForm , corrFun , phiVal , corrForm , subjectID ){
  pgx <- AxioGWAS:::makeLMEPGxPGx2Data( pgx , g[i,],  g2[i,] , NULL , NULL , subVar = subjectID )
  correlationForm <- parse( text = AxioGWAS:::getCorList( list( corrForm = corrForm , phiVal = phiVal , corrFun = corrFun ) ) )
  newForm <- updateModel_2g( pgx , Model , interactionTerm )
  #wrapped with a "try"
  lmTmp <- try(eval( parse( text = paste( "lme( fixed = newForm , data = pgx , random =" , randomForm , ", correlation =" , correlationForm , ", na.action = na.exclude )" ) ) ) )
  
  # if it fails, write why it failed for both regression and anova objects, set pval returned to "NA"
  if(inherits(lmTmp,"try-error")){
    lmTmp<-as.character(lmTmp)
    aTable<-lmTmp
    pval<-NA
  }else{
    aTable <- Anova( lmTmp , type = "III", test.statistic = "Wald")
    lmTmp<-stripLME_callEnv( lmTmp )
    
    if ( is.null( interactionTerm ) ){
      pval<-aTable$"p-value"[which( rownames( aTable ) ==   "Add:Add2"  )]
    }else{
      pval<-aTable$"p-value"[which( rownames( aTable ) ==  paste( interactionTerm , "Add" , sep = ":" ) )]
    }
    
  }
  
  #save values to tables
  writeObjectToTable( lmTmp , paste0( "LME_2Genotypes_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( aTable , paste0( "LME_2Genotypes_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  
  return( list(  pVal = pval , ID=rownames( g )[i] ) )
}



#' run a GLME for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a GLME for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param family The family the generalize model to use. see glm.
#' @param control The control function to use
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a vector of P values, and saves lm models to the db specified in DBPath
#' @export
gwasGLMEFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , family="binomial" , subjectID )
{
  return( lapply( i , gwasGLMEFuncsf , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry , tempFile , family , subjectID ) )
}

#' run a LME for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LMEH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param family The family the generalize model to use. see glm.
#' @param control The control function to use
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P value, and saves lme models to the db specified in dbPath
#'
#' @import lme4
#' @import AxioSerializer
gwasGLMEFuncsf <- function( i , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() , family="binomial" , subjectID )
{
  pgx <- AxioGWAS:::makeLMEPGxData( pgx , g[i,] , NULL , NULL , subVar = subjectID )
  newForm <- updateModel( pgx , Model , interactionTerm )
  
  lmTmp <- try(eval( parse( text = paste( "lme4::glmer( formula = formula(newForm), data = pgx , family =" , family , ", control = lme4::glmerControl(optimizer='bobyqa'))" ) ) ))
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
   
  if(inherits(lmTmp,"try-error")){
    lmTmp<-as.character(lmTmp)
    aTable<-lmTmp
    pval<-NA
  }else{
    aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
    lmTmp<-stripGLME_callEnv( lmTmp )
    if ( is.null( interactionTerm ) ){
      pval<-aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )]
    }else{
      pval<-aTable$"Pr(>Chisq)"[which( rownames( aTable ) ==  paste( interactionTerm , "Add" , sep = ":" ) )]
    }
  }
  
  #save values to tables
  writeObjectToTable( lmTmp , paste0( "LME_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( aTable , paste0( "LME_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  
  return( list(  pVal = pval , ID=rownames( g )[i] ) )
}

#' run a LME for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a LMEH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomForm The name of the random varable to include in results
#' @param correlationForm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P value, and saves lme models to the db specified in dbPath
#'
#' @import nlme
#' @import AxioSerializer
gwasLMEFuncFF <- function( i , g , pgx , Model , interactionTerm , randomForm , correlationForm , subjectID , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- makeLMEPGxData( pgx , g[i,] , NULL , NULL , subVar = subjectID )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- lme( formula( newForm ) , data = pgx , random = randomForm , correlation = correlationForm )
  aTable <- anova( lmTmp )
  writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "LME_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( stripGLM_callEnv( car::Anova( lmTmp , type = "III" , singular.ok = TRUE ) ) , paste0( "LME_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
#    return( list(  pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , data = aTable ) )
    return( list(  pVal = aTable$"p-value"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
#    return( list(  pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , data = aTable ) )
    return( list(  pVal = aTable$"p-value"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
}

#' run a Generalized Linear Model for a GWAS
#'
#' A function called by runGWAS to perform Generalized Linear Model for each snp
#'
#' @param i what row from SNP to call the genotypes from
#' @param g Matrix of "#/#" genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable The name of the table in the DB to save the results to
#' @param snpInfo Details about each snp. Requires at least a column titled 'Chr', for chromosome of snp
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 2 entries: pVal contains a P-Value, lm contains the linear model used
#'
#' @import stats
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,0,0,1,1,1,0,1),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),var4=c(1,2,1,1,2,1,1,2),var5=c(2,2,2,1,1,1,2,1),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasGLMFuncRes<-gwasGLMFunc(i=1,g=exampleSNP, pgx=examplePGx,Model=exampleModel,interactionTerm="var2",genderVar="gender",family="gaussian",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
#'
#' @export
gwasGLMFunc <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , family = "binomial" , snpInfo )
{
  pgx <- makePGxData( pgx , g[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- glm( newForm , data = pgx , family = family )
  writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "GLM_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  aTable <- car::Anova( lmTmp , type = "III" , test = "Wald" , singular.ok = TRUE )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "GLM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) ) + 1] , ID=rownames( g )[i] ) )
  }
}

#' run a GLM for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a GLM for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar the name in pgx for the gender variable
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param subCoef The name of the sub Coeffects varable to include in results
#' @param subEffect The name of the sub Effects varable to include in results
#' @param dbPath The path to database to save results
#' @param Table The name of the table in the DB to save the results to
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a vector of P values, and saves lm models to the table/db specified in DBPath/Table
#'
#' @import stats
#' @export
gwasGLMFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , family = "binomial" )
{
  return( lapply( i , gwasGLMFuncsf , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = overwriteEntry , tempFile = tempFile , family = family ) )
}

#' run a GLM for a GWAS
#'
#' A helper function called by gwasGLMFuncSF in runGWAS to perform a GLM for each snp
#'
#' @param i index of snp of interest in g
#' @param g matrix object of numeric genotype calls, snp by row, sample by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param maf list of snps that pass the maf test
#' @param hwe list of snps that pass the hwe test#'
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 3 entries, a P value, dataframe of results and the name of the target. LM model is also saved to the db specified in dbPath
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasGLMFuncsf <- function( i , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() , family = "binomial" )
{
  pgx <- makePGxData( pgx , g[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- glm( newForm , data = pgx , family = family )
  writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "GLM_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  aTable <- car::Anova( lmTmp , type = "III" , test = "Wald" , singular.ok = TRUE )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "GLM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    return( list(list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) ))
  }
  else
  {
    return( list(list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )] , ID=rownames( g )[i] ) ))
  }
}

#' run a GLM for a GWAS in parallel or series
#'
#' A function called by runGWAS to perform a GLM for each snp, designed for both SnowFall parallel processing and series.
#'
#' @param i the index of the snp in g to run the lm on
#' @param g matrix object of numeric genotype calls, snp by row, sample by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the genderVar in pgx
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable the name of the table to save the lm objects to. done using \code{\link[AxioSerializer]{writeObjectToTable}}
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with two entries, a a P value, and glm model. glm model is also saved to the db specified in databasePath.
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasGLMFuncFF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() , family = "binomial" )
{
  snpVal <- rownames(g)[i]
  
  Add<-g[i,]-1
  if(is.data.frame(Add))
    Add<-do.call('c',Add)
  pgx <- data.frame( pgx , Add = Add )
  
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- glm( newForm , data = pgx , family = family )
  lmWrite <- stripGLM_callEnv( lmTmp )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( lmWrite , paste0( "GLM_" , snpVal ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  aTable <- car::Anova( lmTmp , type = "III" , test = "Wald" , singular.ok = TRUE )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "GLM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    pVal <- aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )]
  }
  else
  {
    iTerm <- c( paste( interactionTerm , "Add" , sep = ":" ) , paste( "Add" , interactionTerm , sep = ":" ) )
    pVal <- aTable$"Pr(>Chisq)"[which( rownames( aTable ) == iTerm[which( iTerm %in% rownames( aTable ) )] )]
  }
  x <- list( pVal = pVal , ID=rownames( g )[i] )
  return( x )
}

#' run a MGLM for a GWAS
#'
#' A function called by runGWAS to perform MGLM for each snp
#'
#' @param i what row from SNP to call the genotypes from
#' @param g Matrix of "#/#" genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable The name of the table in the DB to save the results to
#' @param snpInfo Details about each snp. Requires at least a column titled 'Chr', for chromosome of snp
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 2 entries: pVal contains a P-Value, lm contains the linear model used
#'
#' @import stats
#' @import bit
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,0,0,1,1,1,0,1),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),var4=c(1,2,1,1,2,1,1,2),var5=c(2,2,2,1,1,1,2,1),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasLMEHFuncRes<-gwasMGLMFunc(i=1,g=exampleSNP, pgx=examplePGx,Model=exampleModel,interactionTerm="var2",genderVar="gender",family="gaussian",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
gwasMGLMFunc <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- makePGxData( pgx , g[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- eval( parse( text = paste( "multinom(" , newForm , ", data = pgx , trace = FALSE )" ) ) )
  writeObjectToTable( lmTmp , paste0( "MGLM_" , i ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "MGLM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )] , ID=rownames( g )[i] ) )
  }
}

#' run a MGLM for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a MGLM for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param subCoef The name of the sub Coeffects varable to include in results
#' @param subEffect The name of the sub Effects varable to include in results
#' @param dbPath The path to database to save results
#' @param Table The name of the table in the DB to save the results to
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#'
#' @return a vector of P values, and saves lm models to the db specified in DBPath
#'
#' @import stats
#' @export
gwasMGLMFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  return( sapply( i , gwasMGLMFuncsf , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempFile ) )
}

#' run a MGLM for a GWAS
#'
#' A helper function called by gwasMGLMFuncSF in runGWAS to perform a MGLM for each snp
#'
#' @param i index of snp of interest in g
#' @param g matrix object of numeric genotype calls, snp by row, sample by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param subCoef The name of the sub Coeffects varable to include in results
#' @param subEffect The name of the sub Effects varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 3 entries, a P value, dataframe of results and the name of the target. LM model is also saved to the db specified in dbPath
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasMGLMFuncsf <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- makePGxData( pgx , g[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- eval( parse( text = paste( "multinom(" , newForm , ", data = pgx , trace = FALSE )" ) ) )
  writeObjectToTable( lmTmp , paste0( "MGLM_" , i ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "MGLM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )] , ID=rownames( g )[i] ) )
  }
}

#' run a MGLM for a GWAS in parallel or series
#'
#' A function called by runGWAS to perform a MGLM for each snp, designed for both SnowFall parallel processing and series.
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable the name of the table to save the lm objects to. done using \code{\link[AxioSerializer]{writeObjectToTable}}
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 2 entries, a a P value, and the MGLM model. MGLM model is also saved to the db specified in databasePath.
#'
#' @import stats
#' @import AxioSerializer
#' @import car
#' @export
gwasMGLMFuncFF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- makePGxData( pgx , g[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- eval( parse( text = paste( "multinom(" , newForm , ", data = pgx , trace = FALSE )" ) ) )
  writeObjectToTable( lmTmp , paste0( "MGLM_" , i ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "MGLM_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )] , ID=rownames( g )[i] ) )
  }
}

#' run a coxph for a GWAS
#'
#' A function called by runGWAS to perform coxph for each snp
#'
#' @param i what row from SNP to call the genotypes from
#' @param g Matrix of "#/#" genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable The name of the table in the DB to save the results to
#' @param snpInfo Details about each snp. Requires at least a column titled 'Chr', for chromosome of snp
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 2 entries: pVal contains a P-Value, lm contains the linear model used
#'
#' @import stats
#' @import bit
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,0,0,1,1,1,0,1),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),var4=c(1,2,1,1,2,1,1,2),var5=c(2,2,2,1,1,1,2,1),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' gwasLMEHFuncRes<-gwasMGLMFunc(i=1,g=exampleSNP, pgx=examplePGx,Model=exampleModel,interactionTerm="var2",genderVar="gender",family="gaussian",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' }
#' @export
gwasCOXPHFunc <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- makePGxData( pgx , g$genotypes[i,] , NULL , NULL )
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- coxph( newForm , data = pgx )
  lmWrite <- stripGLM_callEnv( lmTmp )
  aTable <- car::Anova( lmTmp , type = "III" , test.statistic = "Wald" )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( lmWrite , paste0( "COXPH_" , rownames( g$genotypes )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
    writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "COXPH_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )] , ID=rownames( g )[i] ) )
  }
}

#' run a COXPH for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a COXPH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomTerm The name of the random varable to include in results
#' @param heterogeneousTerm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a vector of P values, and saves lm models to the db specified in DBPath
#' @export
gwasCOXPHFuncSF <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() )
{
  return( lapply( i , gwasCOXPHFuncsf , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry , tempFile = tempFile ) )
}

#' run a COXPH for a GWAS in Parallel
#'
#' A function called by runGWAS to perform a COXPH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P value, and saves lme models to the db specified in DBPath
#'
#' @import survival
#' @import AxioSerializer
gwasCOXPHFuncsf <- function( i , g , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  Add<-g[i,]-1
  if(is.data.frame(Add))
    Add<-do.call('c',Add)
  pgx <- data.frame( pgx , Add = Add )
  
  newForm <- updateModel( pgx , Model , interactionTerm )
  lmTmp <- try(coxph( newForm ,data = pgx )) 
  
  # if it fails, write why it failed for both regression and anova objects, set pval returned to "NA"
  if(inherits(lmTmp,"try-error")){
    lmTmp<-as.character(lmTmp)
    aTable<-lmTmp
    pval<-NA
  }else{
    aTable<-coef(summary(lmTmp))
    lmTmp<-stripCOXPH_callEnv( lmTmp )
    if ( is.null( interactionTerm ) ){
      pval<-aTable[which( rownames( aTable ) == "Add" ),"Pr(>|z|)"]
    }else{
      pval<-aTable[which( rownames( aTable ) ==  paste( interactionTerm , "Add" , sep = ":" ) ),"Pr(>|z|)"]
    }
  }
  
  writeObjectToTable(  lmTmp , paste0( "COXPH_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  
  return( list( pVal = pval , ID=rownames( g )[i] ) )
}

#' run a COXPH for a GWAS with two genotypes in Parallel
#'
#' A function called by runGWAS to perform a COXPH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param i what row(s) from SNP to call the genotypes from, input from sflapply(chunk.ffdf(ffobject,by=nrow/ncores),...)
#' @param g ff object of genotype calls, snp by row, subject by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param randomTerm The name of the random varable to include in results
#' @param heterogeneousTerm The name of the heterogeneous varable to include in results
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#' @param g2 object of genotype calls, snp by row, subject by column, for the second genotype to compare

#' @return a vector of P values, and saves lm models to the db specified in DBPath
#' @export
gwasCOXPHFuncSF_2Genotypes <- function( i , g , pgx , Model , interactionTerm , genderVar , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile(), g2 )
{
  return( lapply( i , gwasCOXPHFuncsf_2Genotypes , g , g2 , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry , tempFile = tempFile ) )
}

#' run a COXPH for a GWAS with two genotypes in Parallel
#'
#' A function called by runGWAS to perform a COXPH for each snp, designed for SnowFall parallel processing, and input snp is a  "CNSet" or "SnpSuperSet"
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param g2 object of genotype calls, snp by row, subject by column, for the second genotype to compare
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param dbPath The path to database to save results
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return a P value, and saves lme models to the db specified in DBPath
#'
#' @import survival
#' @import AxioSerializer
gwasCOXPHFuncsf_2Genotypes <- function( i , g , g2 , pgx , Model , interactionTerm , dbPath , dbTable , overwriteEntry = FALSE , tempFile = tempfile() )
{
  Add<-g[i,]-1
  if(is.data.frame(Add))
    Add<-do.call('c',Add)
  Add2<-g2[i,]-1
  if(is.data.frame(Add2))
    Add2<-do.call('c',Add2)
  
  pgx <- data.frame( pgx , Add = Add , Add2 = Add2)
  
  newForm <- updateModel_2g( pgx , Model , interactionTerm )
  lmTmp <- try(coxph( newForm ,data = pgx )) 
  
  # if it fails, write why it failed for both regression and anova objects, set pval returned to "NA"
  if(inherits(lmTmp,"try-error")){
    lmTmp<-as.character(lmTmp)
    aTable<-lmTmp
    pval<-NA
  }else{
    aTable<-coef(summary(lmTmp))
    lmTmp<-stripCOXPH_callEnv( lmTmp )
    if ( is.null( interactionTerm ) ){
      pval<-aTable[which( rownames( aTable ) == "Add:Add2" ),"Pr(>|z|)"]
    }else{
      pval<-aTable[which( rownames( aTable ) ==  paste( interactionTerm , "Add" , sep = ":" ) ),"Pr(>|z|)"]
    }
  }
  
  writeObjectToTable(  lmTmp , paste0( "COXPH_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  
  return( list( pVal = pval , ID=rownames( g )[i] ) )
}

#' run a COXPH for a GWAS in parallel or series
#'
#' A function called by runGWAS to perform a MGLM for each snp, designed for both SnowFall parallel processing and series.
#'
#' @param g vector object of numeric genotype calls for a single snp
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable the name of the table to save the lm objects to. done using \code{\link[AxioSerializer]{writeObjectToTable}}
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 2 entries, a a P value, and the MGLM model. MGLM model is also saved to the db specified in databasePath.
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasCOXPHFuncFF <- function( g , pgx , Model , interactionTerm , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- data.frame( pgx , Add = g[i,] - 1 )
  newForm <- updateModel( pgx , Model , interactionTerm )
#  lmTmp <- coxph( newForm , data = pgx )
  lmTmp <- eval( parse( text = paste( "coxph(" , newForm , ", data = pgx )" ) ) )
  aTable <- car::Anova( lmTmp , type = "III" , test.statistic = "Wald" )
  writeObjectToTable( stripGLM_callEnv( lmTmp ) , paste0( "COXPH_" , rownames( g )[i] ) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "COXPH_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == "Add" )] , ID=rownames( g )[i] ) )
  }
  else
  {
    return( list( pVal = aTable$"Pr(>Chisq)"[which( rownames( aTable ) == paste( interactionTerm , "Add" , sep = ":" ) )] , ID=rownames( g )[i] ) )
  }
}

#' run a Ridge Regression for a GWAS in parallel or series
#'
#' A function called by runGWAS to perform a GLM for each snp, designed for both SnowFall parallel processing and series. Adds interaction term using updateModelRR. Performs RidgeRegression.
#'
#' @param i the index of the snp in g to run the lm on
#' @param g matrix object of numeric genotype calls, snp by row, sample by column
#' @param pgx The data.frame containing the phenotypes used in the model
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the genderVar in pgx
#' @param family The family of Generalized Linear Model to use. Defaults to "binomial", more can be found at \code{\link[stats]{family}}
#' @param dbPath The path to database to save results
#' @param dbTable the name of the table to save the lm objects to. done using \code{\link[AxioSerializer]{writeObjectToTable}}
#' @param overwriteEntry Boolean, passed to \code{\link[AxioSerializer]{writeObjectToTable}} whether to overwrite saved objects in the table specified
#'
#' @return list with 3 entries, a a P value, glm model, and ridge regression coefficients. glm model is also saved to the db specified in databasePath.
#'
#' @import stats
#' @import AxioSerializer
#' @import car
gwasRRFuncFF <- function( i , g , pgx , Model , interactionTerm , genderVar , family = "binomial" , dbPath = NULL , dbTable = NULL , overwriteEntry = FALSE , tempFile = tempfile() )
{
  pgx <- data.frame( pgx , Add = g[i,] - 1 )
  newFormList <- updateModelRR( pgx , Model , interactionTerm )
  pgx <- newFormList$x
  newForm <- newFormList$newForm
  rrPgx <- tryCatch( summary( logisticRidge( newForm , data = pgx , scaling = "none" , lambda = 0.1 ) ) , error = function( e ) return( NULL ) )
  lmWrite <- stripGLM_callEnv( lmTmp )
  if ( !is.null( dbPath ) && !is.null( dbTable ) )
  {
    writeObjectToTable( lmWrite , paste0( "RR_" , rownames( g )[i]) , tableName = dbTable , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  }
  aTable <- car::Anova( lmTmp , type = "III" , singular.ok = TRUE )
  writeObjectToTable( stripGLM_callEnv( aTable ) , paste0( "RR_" , rownames( g )[i] ) , tableName = paste( "Anova" , dbTable , sep = "_" ) , dbPath , overwriteEntry = overwriteEntry , lockFile = tempFile )
  if ( is.null( interactionTerm ) )
  {
    iTerm <- c( paste( interactionTerm , "Add" , sep = ":" ) , "Add" )
  }
  else
  {
    iTerm <- c( paste( interactionTerm , "Add" , sep = ":" ) , paste( "Add" , interactionTerm , sep = ":" ) )
  }
  pVal <- aTable$"Pr(>Chisq)"[which( rownames( aTable ) == iTerm[which( iTerm %in% rownames( aTable ) )] )]
  if ( is.null( rrPgx ) )
  {
    rrPgxCoef <- NULL
  }
  else
  {
    rrPgxCoef <- rrPgx$summaries[["summary1"]]$coefficients[grep( ":" , rownames( rrPgx$summaries[["summary1"]]$coefficients ) , fixed = TRUE),]
  }
  x <- list( pVal = pVal , lm = lmWrite , coef = rrPgxCoef )
  #names( x ) <- rownames( g )[i]
  return( x )
}


