#' Merge Duplicates
#'
#' Merge duplicated values together
#'
#' @param x Sample Dataframe, phenotypes
#' @param g the function to peform the gwas
#' @param uIndx index of unique IDs
#' @param dIndx index of duplicated IDs
#'
#' @return genotype matrix
#'
#' @export
mergeDuplicates <- function( x , g , uIndx = "PAT_ID" , dIndx = "SPCM" )
{
  for ( i in unique( x[,uIndx] ) )
  {
    rNames <- rownames( x[which( x[,uIndx] %in% i ),] )
    genotypes <- g[which( rownames( g ) %in% rNames ),]
    newGenotypes <- matrix( apply( genotypes , 2 , mergeGenotypes ) , nrow = 1 )
    g <- g[-which( rownames( g ) %in% rNames[2] ),]
    g[rNames[1],] <- newGenotypes
  }
  return( g )
}

#' Merge Duplicated Data
#'
#' Merge duplicated values together
#'
#' @param x Sample Dataframe, phenotypes
#' @param XNames the Names of duplicated specimins to merge in dIndx column
#' @param g the function to peform the gwas
#' @param uIndx index of unique IDs
#' @param dIndx index of duplicated IDs
#'
#' @return data.frame
mergeDuplicateData <- function( x , xNames , uIndx = "PAT_ID" , dIndx = "SPCMN" )
{
  xNames <- x[which( x[,dIndx] %in% xNames ),dIndx]
  tMat <- as.data.frame( matrix( nrow = length( xNames ) , ncol = ncol( x ) ) )
  rownames( tMat ) <- xNames
  colnames( tMat ) <- colnames( x )
  for ( i in unique( x[,uIndx] ) )
  {
    tmp <- x[which( x[,uIndx] %in% i ),]
    Indx <- tmp[,dIndx][which( tmp[,dIndx] %in% xNames )]
    tMat[Indx,] <- tmp[Indx,]
    for ( j in names( tmp )[which( !( names( tmp ) %in% c( uIndx , dIndx ) ) )] )
    {
      tMat[Indx,j] <- mergeDuplicateVariable( tmp[,j] )
    }
  }
  return( tMat )
}



#' Merge Duplicated variables
#'
#' Merge duplicated values together
#'
#' @param x  dataframe of variables
mergeDuplicateVariable <- function( x )
{
  if ( inherits( x , c( "factor" , "Date" ) ) )
  {
    x <- as.character( x )
  }
  if ( all( is.na( x ) ) )
  {
    return( NA )
  }
  else if ( setequal( x[1] , x[2] ) )
  {
    return( x[1] )
  }
  else
  {
    idx <- which( !is.na( x ) )
    if ( length( idx ) > 1 )
    {
      cat( "Discordant Data" , x , "\n" )
    }
    else
    {
      return( x[idx] )
    }
  }
}

#' Merge Duplicated Genotuypes
#'
#' Merge duplicated values together
#'
#' @param x dataframe of variables
#'
#' @return data.frame
#'
#' @export
mergeGenotypes <- function( x )
{
  if ( all( is.na( x ) ) )
  {
    return( NA )
  }
  else if ( any( is.na( x ) ) )
  {
    return( x[which( !is.na( x ) )] )
  }
  else
  {
    if ( setequal( x[1] , x[2] ) )
    {
      return( x[1] )
    }
    else
    {
      return( NA )
    }
  }
}


#' Find Discordant Genotypes
#'
#' Merge duplicated values together
#'
#' @param x Sample Dataframe, phenotypes
#' @param g matric of genotypes, snp by row, sample by column
#' @param XNames the Names of duplicated specimins to merge in dIndx column
#' @param g the function to peform the gwas
#' @param uIndx index of unique IDs
#' @param dIndx index of duplicated IDs
#'
#' @return data.frame
#'
#' @export
findDiscordantCalls <- function( x , g , uIndx = "PAT_ID" , dIndx = "SPCM" )
{
  discords <- matrix( NA , nrow = length( unique( x[,uIndx] ) ) , ncol = ncol( g ) )
  rownames( discords ) <- unique( x[,uIndx] )
  for ( i in unique( x[,uIndx] ) )
  {
    rNames <- rownames( x[which( x[,uIndx] %in% i ),] )
    genotypes <- g[which(colnames( g ) %in% rNames ),]
    discords[as.character( i ),] <- matrix( apply( genotypes , 2 , getDiscordantCall ) , nrow = 1 )
  }
  return( colSums( discords ) )
}

#' Get Discordant Genotypes
#'
#' determine if calls are discordant or not
#'
#' @param x vector with 2 values
#'
#' @return data.frame
#'
#' @export
getDiscordantCall <- function( x )
{
  if ( all( is.na( x ) ) )
  {
    return( NA )
  }
  else if ( any( is.na( x ) ) )
  {
    return( TRUE )
  }
  else
  {
    if ( setequal( x[1] , x[2] ) )
    {
      return( FALSE )
    }
    else
    {
      return( TRUE )
    }
  }
}
#' get IBS
#'
#' calculate identical by state value
#'
#' @param g matrix of genotypes, snp by row, sample by column
#' @param snowFall Boolean value for whether to perform this analysis in parallel or not. defaults to FALSE
#'
#' @return matrix of IBS values
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=9,byrow=TRUE)
#'  rownames(exampleSNP)<-c(1,2,3,4,5,6,7,8,9)
#'  
#'  getIBS(exampleSNP,FALSE)
#' }
#' @export
getIBS <- function( g , snowFall = FALSE )
{
  ibsMat <- diag( nrow( g ) )
  pNames <- rownames( g )
  rownames( ibsMat ) <- pNames
  colnames( ibsMat ) <- pNames
  while ( length( pNames ) > 1 )
  {
    pName <- pNames[1]
    pNames <- setdiff( pNames , pName )
    indx1 <- which( rownames( g ) %in% pName )
    indx2 <- which( rownames( g ) %in% pNames )
    if ( snowFall )
    {
      if ( length( indx2 ) > 1 )
      {
        ibsMat[indx1,indx2] <- ibsMat[indx2,indx1] <- unlist( sfApply( cbind( rep( indx1 , length( indx2 ) ) , indx2 ) , 1 , function( x , g ) { sum( g[x[1],] == g[x[2],] , na.rm =TRUE ) / sum( unlist( apply( g[x,] , 2 , function( x ) any( !is.na( x ) ) ) ) ) } , g ) )
      }
      else
      {
        ibsMat[indx1,indx2] <- ibsMat[indx2,indx1] <- unlist( apply( cbind( rep( indx1 , length( indx2 ) ) , indx2 ) , 1 , function( x , g ) { sum( g[x[1],] == g[x[2],] , na.rm =TRUE ) / sum( unlist( apply( g[x,] , 2 , function( x ) any( !is.na( x ) ) ) ) ) } , g ) )
      }
    }
    else
    {
      ibsMat[indx1,indx2] <- ibsMat[indx2,indx1] <- unlist( apply( cbind( rep( indx1 , length( indx2 ) ) , indx2 ) , 1 , function( x , g ) { sum( g[x[1],] == g[x[2],] , na.rm =TRUE ) / sum( unlist( apply( g[x,] , 2 , function( x ) any( !is.na( x ) ) ) ) ) } , g ) )
    }
  }
  return( ibsMat )
}

#' get IBS
#'
#' calculate identical by state value
#'
#' @param x matrix of IBS values from \code{getIBS}
#' @param cutOff Level at which to identify subjects with IBS > cutoff
#'
#' @return matrix of IBS values
#' 
#' @examples
#'  exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=6,byrow=TRUE)
#'  rownames(exampleSNP)<-c(1,2,3,4,5,6)
#'  
#'  ibsMat<-getIBS(exampleSNP,FALSE)
#'
#'  summaryIBSMat(ibsMat,.75)
#'  
#' @export
summaryIBSMat <- function( x , cutOff = 0.95 )
{
  indx <- which( ( x - diag( nrow( x ) ) ) > cutOff , arr.ind = TRUE )
  cat( nrow( indx ) , "individuals have an IBS coefficient >" , cutOff , "\n" )
  cat( "The range of IBS coefficients is: Low" , min( x ) , ", Max" , max( x - diag( nrow( x ) ) ) , "\n" )
  return( indx )
}
#' calculate IBS
#'
#' calculates IBS number for two calls, either 1 or 0
#'
#' @param x vector containing 2 values
#'
#' @export
getIBSNum <- function( x )
{
  return( as.integer( setequal( x[1] , x[2] ) ) )
}

#' calculate IBS
#'
#' calculates IBS denominator
#'
#' @param x vector
#'
#' @export
getIBSDenom <- function( x )
{
  return( as.integer( any( !is.na( x ) ) ) )
}

#----------------------------------------------------------------

loadSNPMatrix <- function( x )
{
  rNames <- rownames( x )
  cNames <- colnames( x )
  g <- as.data.frame( matrix( unlist( lapply( x , function( x ) makeSNPMatrix( x ) ) ) , nrow = nrow( x ) ) )
  g <- matrix( unlist( lapply( g , function( x ) as.numeric( x ) ) ) , ncol = ncol( g ) )
  g[which( is.na( g ) , arr.ind = TRUE )] <- 0
  colnames( g ) <- cNames
  rownames( g ) <- rNames
  print( g[1:5,1:10] )
  g <- new( "SnpMatrix" , g )
  return( g )
}

makeSNPMatrix <- function( x )
{
  if ( all( is.na( x ) ) )
  {
    return( x )
  }
  else
  {
    alleles <- table( unlist( lapply( as.character( x ) , function( x ) strsplit( x , "/" , fixed = TRUE ) ) ) )
    if ( length( alleles ) < 2 )
    {
      factLevels <- c( paste( names( alleles ) , names( alleles ) , sep = "/" ) )
    }
    else
    {
      if ( isTRUE( all.equal( as.integer( alleles[1] ) , as.integer( alleles[2] ) ) ) )
      {
        if ( length( table( x ) ) > 2 )
        {
          factLevels <- names( table( x ) )
        }
        else
        {
          homwt <- paste( names( alleles )[1] , names( alleles )[1] , sep = "/" )
          het <- paste( names( alleles )[1] , names( alleles )[2] , sep = "/" )
          hommt <- paste( names( alleles )[2] , names( alleles )[2] , sep = "/" )
          factLevels <- c( homwt , het , hommt )
        }
      }
      else
      {
        homwt <- paste( names( alleles )[which.max( alleles )] , names( alleles )[which.max( alleles )] , sep = "/" )
        het <- paste( names( alleles )[which.max( alleles )] , names( alleles )[which.min( alleles )] , sep = "/" )
        hommt <- paste( names( alleles )[which.min( alleles )] , names( alleles )[which.min( alleles )] , sep = "/" )
        factLevels <- c( homwt , het , hommt )
      }
    }
    return( factor( x , levels = factLevels ) )
  }
}

updateCoefficients <- function( x , pgx , subCoef , subEffect , rsid , ffdf = FALSE , Model = NULL )
{
  vars <- attr( terms( formula( Model ) ) , "term.labels" )
  pgx <- data.frame( pgx , Add = runif( nrow( pgx ) , 0 , 2 ) )
  cfnts <- makeCoefLabels( vars , pgx )
  dfnames <- c( "Raw P-Value" , paste( "Estimate" , cfnts ) , paste( "Std. Err." , cfnts ) , paste( "p-value" , cfnts ) , "Formula" )
  dataFrame <- resultsDataFrame( Model , dfnames )
  rownames( dataFrame ) <- rsid
  if ( !is.null( subCoef ) )
  {
    indx <- c( 1 )
    for ( i in subCoef )
    {
      indx <- c( indx , grep( i , names( dataFrame ) , fixed = TRUE ) )
    }
    dataFrame <- dataFrame[,indx]
  }
  if ( !is.null( subEffect ) )
  {
    indx <- c( 1 )
    for ( i in subEffect )
    {
      indx <- c( indx , grep( i , names( dataFrame ) , fixed = TRUE ) )
    }
    dataFrame <- dataFrame[,indx]
  }
  if ( is.null( x ) )
  {
    x <- as.ffdf( dataFrame )
    finalizer( x ) <- "close"
  }
  else
  {
    x <- update.ffdf( x , dataFrame )
  }
  return( x )
}

getCoefficients <- function( x , pgx , subCoef , subEffect , ffdf = FALSE )
{
  if ( length( x$modelFits ) > 1 )
  {
    vars <- attr( terms( formula( x$modelFits[[1]]$lm ) ) , "term.labels" )
    depVar <- rownames( attr( terms( formula( x$modelFits[[1]]$lm ) ) , "factors" ) )[attr( terms( formula( x$modelFits[[1]]$lm ) ) , "response" )]
  }
  else
  {
    vars <- attr( terms( formula( x$modelFits[[names( x$modelFits )]]$lm ) ) , "term.labels" )
    depVar <- rownames( attr( terms( formula( x$modelFits[[names( x$modelFits )]]$lm ) ) , "factors" ) )[attr( terms( formula( x$modelFits[[names( x$modelFits )]]$lm ) ) , "response" )]
  }
  pgx <- data.frame( pgx , Add = runif( nrow( pgx ) , 0 , 2 ) )
  if ( inherits( pgx[,depVar] , "factor" ) )
  {
    if ( length( levels( pgx[,depVar] ) ) > 2 )
    {
      levs <- levels( pgx[,depVar] )[-1]
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
  cfnts <- makeCoefLabels( vars , pgx , levs )
  dfnames <- c( "Raw P-Value" , paste( "Estimate" , cfnts ) , paste( "Std. Err." , cfnts ) , paste( "p-value" , cfnts ) , "Formula" )
  if ( ffdf )
  {
    dataFrame <- data.frame( x$adjPVals , matrix( resultsDataFrame( x$modelFits[[names( x$modelFits )]] , dfnames ) , nrow = 1 ) , stringsAsFactors = FALSE )
  }
  else
  {
    tmp <- lapply( x$modelFits , function( x , dfnames ) { resultsDataFrame( x , dfnames ) } , dfnames )
    dataFrame <- data.frame( x$adjPVals , matrix( unlist( tmp ) , nrow = length( x$adjPVals ) , byrow = TRUE ) )
  }
  rownames( dataFrame ) <- names( x$modelFits )
  names( dataFrame ) <- c( "Adjusted P-Value" , dfnames )
  if ( !is.null( subCoef ) )
  {
    indx <- c( 1 , 2 )
    for ( i in subCoef )
    {
      indx <- c( indx , grep( i , names( dataFrame ) , fixed = TRUE ) )
    }
    dataFrame <- dataFrame[,indx]
  }
  if ( !is.null( subEffect ) )
  {
    indx <- c( 1 , 2 )
    for ( i in subEffect )
    {
      indx <- c( indx , grep( i , names( dataFrame ) , fixed = TRUE ) )
    }
    dataFrame <- dataFrame[,indx]
  }
  return( dataFrame )
}

getCoefficients2 <- function( x , pgx , subCoef , subEffect , ffdf = FALSE )
{
  lmTmp <- x$modelFits[[1]][[1]]$lm
  vars <- attr( terms( formula( lmTmp ) ) , "term.labels" )
  depVar <- rownames( attr( terms( formula( lmTmp ) ) , "factors" ) )[attr( terms( formula( lmTmp ) ) , "response" )]
  if ( !is.null( x$modelFits[[1]][[1]]$coef ) )
  {
    depVar <- gsub( "Num" , "" , depVar , fixed = TRUE )
  }
  pgx <- data.frame( pgx , Add = runif( nrow( pgx ) , 0 , 2 ) )
  if ( inherits( pgx[,depVar] , "factor" ) )
  {
    if ( length( levels( pgx[,depVar] ) ) > 2 )
    {
      levs <- levels( pgx[,depVar] )[-1]
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
  cfnts <- makeCoefLabels( vars , pgx , levs )
  dfnames <- c( "Raw P-Value" , paste( "Estimate" , cfnts ) , paste( "Std. Err." , cfnts ) , paste( "p-value" , cfnts ) , "Formula" )
  tmp <- lapply( x$modelFits , function( x , dfnames ) { resultsDataFrame( x[[1]] , dfnames ) } , dfnames )
  dataFrame <- data.frame( x$adjPVals , matrix( unlist( tmp ) , nrow = length( x$adjPVals ) , byrow = TRUE ) )
  rownames( dataFrame ) <- unlist( lapply( x$modelFits , function( x ) names( x ) ) )
  names( dataFrame ) <- c( "Adjusted P-Value" , dfnames )
  if ( !is.null( subCoef ) )
  {
    indx <- c( 1 , 2 )
    for ( i in subCoef )
    {
      indx <- c( indx , grep( i , names( dataFrame ) , fixed = TRUE ) )
    }
    dataFrame <- dataFrame[,indx]
  }
  if ( !is.null( subEffect ) )
  {
    indx <- c( 1 , 2 )
    for ( i in subEffect )
    {
      indx <- c( indx , grep( i , names( dataFrame ) , fixed = TRUE ) )
    }
    dataFrame <- dataFrame[,indx]
  }
  return( dataFrame )
}

getRRCoefficients <- function( x )
{
  lmTmp <- x$modelFits[[1]][[1]]$coef
  cfnts <- rownames( lmTmp )
  dfnames <- c( paste( "Ridge Estimate" , cfnts ) , paste( "Ridge Std. Err." , cfnts ) , paste( "Ridge p-value" , cfnts ) )
  tmp <- lapply( x$modelFits , function( x ) { cbind( x[[1]]$coef[,"Estimate"] , x[[1]]$coef[,"Std. Error"] , x[[1]]$coef[,"Pr(>|t|)"] ) } )
  dataFrame <- data.frame( matrix( unlist( tmp ) , nrow = length( x$adjPVals ) , byrow = TRUE ) )
  rownames( dataFrame ) <- unlist( lapply( x$modelFits , function( x ) names( x ) ) )
  names( dataFrame ) <- dfnames
  return( dataFrame )
}

writeResults <- function( x , g , file , pgx , sheet , sort = NULL , subCoef = NULL , subEffect = NULL , MAF = FALSE , HWE = FALSE , hweAdjust = "none" )
{
  ffdf <- inherits( g , c( "CNSet" , "SnpSuperSet" ) )
  dataFrame <- getCoefficients( x , pgx , subCoef , subEffect , ffdf )
  if ( MAF )
  {
    cNames <- names( dataFrame )
    maf <- unlist( apply( g$genotypes , 2 , mafFreq ) )
    dataFrame <- data.frame( maf = maf[rownames( dataFrame )] , dataFrame )
    names( dataFrame ) <- c( "MAF" , cNames )
  }
  if ( HWE )
  {
    cNames <- names( dataFrame )
    hwe <- hweTest( g , adjust = hweAdjust )
    dataFrame <- data.frame( hwe = hwe[rownames( dataFrame )] , dataFrame )
    names( dataFrame ) <- c( "HWE" , cNames )
  }
  if ( !is.null( g$snpInfo ) )
  {
    cNames <- names( dataFrame )
    dataFrame <- data.frame( chromosome = g$snpInfo[rownames( dataFrame ),c( "Chr" , "GeneName" )] , dataFrame )
    names( dataFrame ) <- c( "Chromosome" , "Gene Name" , cNames )
  }
  if ( !is.null( sort ) )
  {
    dataFrame <- switch( sort ,
        Raw = dataFrame[order( dataFrame[,"Raw P-Value"] ),] ,
        Adjusted = dataFrame[order( dataFrame[,"Adjusted P-Value"] ),]
    )
  }
  WB <- loadWorkbook( file.path( file ) , create = TRUE )
  createSheet( WB , name = sheet )
  writeWorksheet( WB , dataFrame , sheet = sheet , header = TRUE , rownames = "rownames" )
  saveWorkbook( WB )
}

writeResults2 <- function( x , g , gs , file , pgx , sheet , sort = NULL , subCoef = NULL , subEffect = NULL , MAF = FALSE , HWE = FALSE , hweAdjust = "none" )
{
  ffdf <- inherits( g , c( "CNSet" , "SnpSuperSet" ) )
  dataFrame <- getCoefficients2( x , pgx , subCoef , subEffect , ffdf )
  if ( MAF )
  {
    cNames <- names( dataFrame )
    maf <- unlist( apply( calls( g )[] , 1 , mafFreqIllumina , chromosome( g ) , pData( g )[,"Sex"] ) )
    dataFrame <- data.frame( maf = maf[rownames( dataFrame )] , dataFrame )
    names( dataFrame ) <- c( "MAF" , cNames )
  }
  if ( HWE )
  {
    cNames <- names( dataFrame )
    hwe <- hweTest( g , adjust = hweAdjust )
    dataFrame <- data.frame( hwe = hwe[rownames( dataFrame )] , dataFrame )
    names( dataFrame ) <- c( "HWE" , cNames )
  }
  if ( ffdf )
  {
    cNames <- names( dataFrame )
    chr <- chromosome( g )
    names( chr ) <- featureNames( g )
    dataFrame <- data.frame( chromosome = chr[rownames( dataFrame )] , dataFrame )
    names( dataFrame ) <- c( "Chromosome" , cNames )
  }
  if ( !is.null( g$snpInfo ) )
  {
    cNames <- names( dataFrame )
    dataFrame <- data.frame( chromosome = g$snpInfo[rownames( dataFrame ),c( "Chr" , "GeneName" )] , dataFrame )
    names( dataFrame ) <- c( "Chromosome" , "Gene Name" , cNames )
  }
  if ( !is.null( sort ) )
  {
    dataFrame <- switch( sort ,
        Raw = dataFrame[order( dataFrame[,"Raw P-Value"] ),] ,
        Adjusted = dataFrame[order( dataFrame[,"Adjusted P-Value"] ),]
    )
  }
  if ( !is.null( x$modelFits[[1]][[1]]$coef ) )
  {
    cNames <- names( dataFrame )
    dataTmp <- getRRCoefficients( x )
    dataFrame <- data.frame( dataFrame , dataTmp[rownames( dataFrame ),] )
    names( dataFrame ) <- c( cNames , names( dataTmp ) )
  }
  xNames <- names( dataFrame )
  dataFrame <- data.frame( "Gene Name" = gs[rownames( dataFrame )] , dataFrame )
  names( dataFrame ) <- c( "Gene Name" , xNames )
  WB <- loadWorkbook( file.path( file ) , create = TRUE )
  createSheet( WB , name = sheet )
  writeWorksheet( WB , dataFrame , sheet = sheet , header = TRUE , rownames = "rownames" )
  saveWorkbook( WB )
}

writeResultsSNPWise <- function( x , g , pgx , MAF = NULL , HWE = NULL , subCoef = NULL , subEffect = NULL )
{
  dataFrame <- getCoefficients( x , pgx , subCoef , subEffect , TRUE )
  if ( !is.null( MAF ) )
  {
    cNames <- names( dataFrame )
    dataFrame <- data.frame( MAF , dataFrame )
    names( dataFrame ) <- c( "MAF" , cNames )
  }
  if ( !is.null( HWE ) )
  {
    cNames <- names( dataFrame )
    dataFrame <- data.frame( HWE , dataFrame )
    names( dataFrame ) <- c( "HWE" , cNames )
  }
  if ( !is.null( g$snpInfo ) )
  {
    cNames <- names( dataFrame )
    dataFrame <- data.frame( chromosome = g$snpInfo[rownames( dataFrame ),c( "Chr" , "GeneName" )] , dataFrame )
    names( dataFrame ) <- c( "Chromosome" , "Gene Name" , cNames )
  }
  return( dataFrame )
}

resultsDataFrame <- function( x , dfnames )
{
  if ( inherits( x$lm , "multinom" ) )
  {
    tmp <- summary( x$lm )
    cfnts <- cbind( as.vector( tmp$coefficients ) , as.vector( tmp$standard.errors ) )
    if ( length( tmp$lev ) > 2 )
    {
      rnames <- vector()
      for ( i in rownames( tmp$coefficients ) )
      {
        rnames <- c( rnames , paste( i , tmp$coefnames , sep = "." ) )
      }
      rownames( cfnts ) <- rnames
    }
    else
    {
      rownames( cfnts ) <- tmp$coefnames
    }
    cfnts <- cbind( cfnts , cfnts[,1] / cfnts[,2] )
    cfnts <- cbind( cfnts , pt( cfnts[,3] , 1 ) )
  }
  else
  {
    cfnts <- coefficients( summary( x$lm ) )
  }
  df <- c( x[["pVal"]] , cfnts[,1] , cfnts[,2] , cfnts[,4] , deparse( formula( x$lm ) ) )
  names( df ) <- c( "Raw P-Value" , paste( "Estimate" , rownames( cfnts ) ) , paste( "Std. Err." , rownames( cfnts ) ) , paste( "p-value" , rownames( cfnts ) ) , "Formula" )
  if ( length( df ) < length( dfnames) )
  {
    dfTmp <- df
    df <- rep( NA , length( dfnames ) )
    names( df ) <- dfnames
    df[names( dfTmp )] <- dfTmp
  }
  return( df )
}

resultsDataFrame2 <- function( x , dfnames )
{
  if ( inherits( x$lm , "multinom" ) )
  {
    tmp <- summary( x$lm )
    cfnts <- cbind( as.vector( tmp$coefficients ) , as.vector( tmp$standard.errors ) )
    if ( length( tmp$lev ) > 2 )
    {
      rnames <- vector()
      for ( i in rownames( tmp$coefficients ) )
      {
        rnames <- c( rnames , paste( i , tmp$coefnames , sep = "." ) )
      }
      rownames( cfnts ) <- rnames
    }
    else
    {
      rownames( cfnts ) <- tmp$coefnames
    }
    cfnts <- cbind( cfnts , cfnts[,1] / cfnts[,2] )
    cfnts <- cbind( cfnts , pt( cfnts[,3] , 1 ) )
  }
  else
  {
    cfnts <- coefficients( summary( x$lm ) )
  }
  df <- c( x[["pVal"]] , cfnts[,1] , cfnts[,2] , cfnts[,4] , deparse( formula( x$lm ) ) )
  names( df ) <- c( "Raw P-Value" , paste( "Estimate" , rownames( cfnts ) ) , paste( "Std. Err." , rownames( cfnts ) ) , paste( "p-value" , rownames( cfnts ) ) , "Formula" )
  if ( length( df ) < length( dfnames) )
  {
    dfTmp <- df
    df <- rep( NA , length( dfnames ) )
    names( df ) <- dfnames
    df[names( dfTmp )] <- dfTmp
  }
  return( df )
}

makeCoefLabels <- function( x , pgx , levs = NULL )
{
  newLabs <- "(Intercept)"
  for ( i in x )
  {
    if ( length( grep( ":" , i , fixed = TRUE ) ) > 0 )
    {
      j <- unlist( strsplit( i , ":" , fixed = TRUE ) )
      tvar <- vector( mode = "list" , length = length( j ) )
      names( tvar ) <- j
      for ( k in j )
      {
        if ( inherits( pgx[,k] , "factor" ) )
        {
          tvar[[k]] <- paste( k , levels( pgx[,k] )[-1] , sep = "" )
        }
        else
        {
          tvar[[k]] <- k
        }
      }
      newLabs <- c( newLabs , apply( expand.grid( tvar ) , 1 , function( x ) paste( x , collapse = ":" ) ) )
    }
    else
    {
      if ( inherits( pgx[,i] , "factor" ) )
      {
        newLabs <- c( newLabs , paste( i , levels( pgx[,i] )[-1] , sep = "" ) )
      }
      else
      {
        newLabs <- c( newLabs , i )
      }
    }
  }
  if ( !is.null( levs ) )
  {
    tmp <- vector()
    for ( i in levs )
    {
      tmp <- c( tmp , paste( i , newLabs , sep = "." ) )
    }
    newLabs <- tmp
  }
  return( newLabs )
}

getCorList <- function( x )
{
  if ( inherits( x , "corStruct" ) )
  {
    xList <- strsplit( deparse( summary( x ) ) , split = "," )[[1]]
    phiVal <- tanh( as.numeric( gsub( "structure[[:punct:]]" , "" , xList[1] ) ) / 2 )
    formula <- xList[2]
    formula <- strsplit( formula , split = "=" )[[1]][2]
    fun <- strsplit( xList[6] , split = "\"" )[[1]][2]
    return( list( corrFun = fun , phiVal = phiVal , corrForm = formula ) )
  }
  else
  {
    return( paste( x$corrFun , "( " , x$phiVal , " , form = " , x$corrForm , " )" , sep = "" ) )
  }
}

getCor <- function( x )
{
  if ( inherits( x , "corStruct" ) )
  {
    xList <- strsplit( deparse( summary( x ) ) , split = "," )[[1]]
    phiVal <- tanh( as.numeric( gsub( "structure[[:punct:]]" , "" , xList[1] ) ) / 2 )
    formula <- xList[2]
    formula <- gsub( "formula" , "form" , formula )
    fun <- strsplit( xList[6] , split = "\"" )[[1]][2]
    return( paste( fun , "! " , phiVal , "," , formula , " !" , sep = "" ) )
  }
  else
  {
    return( gsub( "! " , "( " , gsub( " !" , " )" , x ) ) )
  }
}
