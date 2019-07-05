#' Run a contrast on the results of the GWAS
#'
#' Input one of the gwas functions, and it will perform contrasts and return a list with entries for estOut and pValue
#'
#' @param i index of of snp in fd
#' @param func the function to peform the gwas
#' @param pgx DF of phenotypes for each subject
#' @param fd Matrix of genotype calls, snp by row, subject by column. Genotype format varies by called "func".
#' @param Model the original formula to update. Must be a formula or model.frame (required)
#' @param interactionTerm One of the colnames of x, to be added to the model, and interact with the genotype
#' @param genderVar The name of the gender varable in pgx
#' @param con contrasts to perform. You can perform multiple contrasts in one. Takes the form: list(phenotype=list(pos=c("phenotypeLevel"),neg=c("phenotypeLevel")),...)
#' @param ... additional arguments to be passed to func
#'
#' @return a list object
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#' exampleCon<-list( "gender:var2" = list(pos="Male/2",neg=c("Male/3","Female/2","Female/3")))
#'
#' newContrast<-runContrastGWAS(1,gwasFunc,exampleSNP,examplePGx,exampleModel,"var2","gender",exampleCon,
#'                databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#'}
runContrastGWAS <- function( i , func , fd , pgx , Model , interactionTerm , genderVar , con , ... )
{
  lmOut <- func( i , fd , pgx , Model , interactionTerm , genderVar , ... )$lm
  return( contrastGWAS( lmOut , con ) )
}

#' Run a contrast on the lmObject
#'
#' Input an lm object and perform a contrast based on your input
#'
#' @param lmOut LM object
#' @param con contrasts to perform. You can perform multiple contrasts in one. Takes the form: list(phenotype=list(pos=c("phenotypeLevel"),neg=c("phenotypeLevel")),...)
#'
#' @return a list object
#'
#' @import doBy
#' 
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#' exampleCon<-list( "gender" = list(pos="Male",neg=c("Female")))
#' 
#' lmObject<-gwasFunc(1,exampleSNP,examplePGx,exampleModel,"var2","gender",databasePath="runGWAS.DB", dbTable="gwasFunc", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#'
#' contrastGWAS(lmObject$lm,exampleCon)
#'}
contrastGWAS <- function( lmOut , con )
{
  lambdaMat <- matrix( 0 , nrow = length( con ) , ncol = length( coef( lmOut ) ) )
  coefNames <- names( coef( lmOut ) )
  naCoefNames <- names( coef( lmOut )[which( is.na( coef( lmOut ) ) )] )
  r <- 1
  for ( i in seq( 1 , length( con ) ) )
  {
    lengthPos <- length( con[[i]]$pos ) - length( which( naCoefNames == paste( names( con[i] ) , con[[i]]$pos , sep = "" ) ) )
    lengthNeg <- length( con[[i]]$neg ) - length( which( naCoefNames == paste( names( con[i] ) , con[[i]]$neg , sep = "" ) ) )
    for ( j in con[[i]]$pos )
    {
      lambdaMat[r,which( coefNames == paste( names( con[i] ) , j , sep = "" ) )] <- ifelse( lengthPos > 0 , 1 / lengthPos , 0 )
    }
    for ( j in con[[i]]$neg )
    {
      lambdaMat[r,which( coefNames == paste( names( con[i] ) , j , sep = "" ) )] <- ifelse( lengthPos > 0 , -1 / lengthNeg , 0 )
    }
    r <- r + 1
  }
  if ( length( naCoefNames ) > 0 )
  {
    lambdaMat <- lambdaMat[,-which( coefNames %in% naCoefNames )]
  }
  estOut <- esticon( lmOut , lambdaMat )
  print( estOut )
  #print( table( lmOut$data[,"Add"] , lmOut$data[,"TreatmentRed"] , lmOut$data[,"ACR20"] ) )
  print( car::vif( lmOut ) )
  return( list( estOut = estOut , pval = estOut[,grep( "Pr(>" , colnames( estOut ) , fixed = TRUE )] ) )
}


#' Run a GWAS and Contrast
#'
#' Perform both a GWAS and Contrast at the same time. Can be performed in parallel(snowfall) or series.
#'
#' @param pgx DF of phenotypes for each subject
#' @param snp Matrix of genotype calls, snp by row, subject by column. Genotype format varies by called "func".
#' @param Model The original base model for the GWAS, excluding snp and snp interaction terms
#' @param func The name of the function to apply in the GWAS
#' @param interactionTerm The name of the variable to be included as the interaction term with the snp
#' @param genderVar The name of the gender varable in pgx
#' @param adjust The method to adjust the P-values by after the GWAS. defaults to "BY", see \code{\link[stats]{p.adjust}} for more options
#' @param replace The name of the variable to remove from the Model if necessary
#' @param snowFall Boolean Value, Should the GWAS be run in parallel or series. If TRUE, use \code{\link[snowfall]{sfInit}} prior to \code{runGWAS}
#' @param file Name of file which to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param subCoeff coefficients to keep and to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param subEffect Effect to keep and to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param conFunc Function to perform contrast. Defaults to \code{runContrastGWAS}.
#' @param con contrasts to perform. You can perform multiple contrasts in one. Takes the form: list(phenotype=list(pos=c("phenotypeLevel"),neg=c("phenotypeLevel")),...)
#' @param dbPath path to database to save results in if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param dbConnection a sqlite connection to a DB, alternative to dbPath to save results if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param dbTABLENAME The name of the table in the DB to save the results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE
#' @param ... additional arguments passed to func 
#'
#' @return a list object
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
#' exampleSNPInfo<-data.frame(Chr=c(1,2))
#' exampleModel<-formula("var1~var2+var3+gender")
#' exampleCon<-list( "gender" = list(pos="Male",neg=c("Female")))
#' 
#' output<-runContrasts(pgx=examplePGx,snp=exampleSNP,Model=exampleModel,func = gwasFunc,interactionTerm="var2",genderVar="gender",databasePath="runGWAS.DB", dbTable="runContrast_Test", snpInfo=exampleSNPInfo, overwriteEntry=TRUE, con = exampleCon)
#'}
runContrasts <- function( pgx , snp , Model , func = gwasFunc , interactionTerm = "ARMACD" , genderVar = "SEX.Decode" , adjust = "BY" , replace = NULL , snowFall = FALSE , file = NULL , subCoef = NULL , subEffect = NULL , conFunc = runContrastGWAS , con = NULL , dbPath = NULL , dbConnection = NULL , dbTABLENAME = NULL , ... )
{
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
      sfExport( "fd" )
      sfClusterEval( open( fd ) )
      sfExport( "gwasFuncsf" )
      sfExport( "makePGxData" )
      sfExport( "makeAdditive" )
      sfExport( "additiveGenotype" )
      sfExport( "updateModel" )
      sfExport( "reducedModel" )
      sfExport( "getCoefficients" )
      sfExport( "makeCoefLabels" )
      sfExport( "resultsDataFrame" )
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
      }
      
      #Create a table with all the columns in dfTmp by Adjusted P-Value to allow merging later
      databaseQuery<-paste("CREATE TABLE",dbTABLENAME,"(SNP varchar ,",paste0('"',setdiff(names(dfTmp),"Adjusted P-Value"),'"',collapse=" double , "),");")
      dbSendQuery(dbConnection,databaseQuery)
      dbDisconnect(dbConnection)
      
      #export function to be sfLappy-ed to each worker
      #eval(parse(text=paste("sfExport('",func,"')")))
      
      ###########
      ## RUN PARALLELIZED LMFIT FUNCTION
      ##########
      gwasRes <- sfLapply( chunk.ffdf( fd , by = nrow( fd ) / length( sfGetCluster() ) ) , conFunc , func , fd , pgx , Model , interactionTerm , subCoef , subEffect , dbPath , dbTABLENAME , con=con, ... )
      
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
      
      write.csv( file = file , modelFits )
      
      return( list( adjPVals = adjustedPValue , modelFits = modelFits ) )
    }
    else
    {
      gwasRes <- lapply( seq( 1 , dim( snp )["Features"] ) , conFunc , func , fd , pgx , Model , interactionTerm , genderVar , con=con , ... )
      close( fd )
      adjPVals <-  p.adjust( unlist( lapply( gwasRes , function( x ) x$pval ) ) , adjust )
      return( list( adjPVals = adjPVals , modelFits = gwasRes ) )
    }
  }
  else
  {
    if ( snowFall )
    {
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "AxioSerializer" , ( .packages() ) , fixed = TRUE ) ) ) ) ) ) )
      {
        sfLibrary( AxioSerializer , keep.source = FALSE )
      }
      
      #added confunc to arguments
      gwasRes <- sfLapply( seq( 1 , nrow( snp ) ), conFunc , func , snp , pgx , Model , interactionTerm , genderVar , con = con , ... )
    }
    else
    {
      #added confunc
      gwasRes <- lapply( seq( 1 , nrow( snp ) ) , conFunc, func , snp , pgx , Model , interactionTerm , genderVar , con = con , ... )
    }
    # names( gwasRes ) <- colnames( snp )
    # adjPVals <-  p.adjust( unlist( lapply( gwasRes , function( x ) x$pVal ) ) , adjust )
    
    adjPVals <-  data.frame(do.call('rbind', lapply( gwasRes , function( x ) c('pval'=as.numeric(as.character(x$pval)) ) )))
    print(adjPVals)
    adjPVals$adjpval <- p.adjust( as.numeric(as.character(adjPVals$pval)) , adjust )
    adjPVals$snp<-rownames(snp)[as.numeric(rownames(adjPVals))]
    
    return( list( adjPVals = adjPVals , modelFits = gwasRes ) )
  }
}
