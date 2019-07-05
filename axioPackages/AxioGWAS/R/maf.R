
#' Calculate MAF
#'
#' Calculate MAF for snp
#'
#' @param i index of Matrix of X
#' @param x Matrix of genotypes
#' @param Chr vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp

#'
#' @return numeric Minor Allele Frequency
#'
#' @examples
#' \dontrun{
#'   exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=4,byrow=TRUE)
#'
#'  mafFreq(1,exampleSNP,c(1,1,2,2),c("M","M","F","M","F","F","F","F","M"))
#' }
#' @export
mafFreq <- function( i , x , Chr = NULL , gender = NULL )
{
  if ( all( is.na( x[i,] ) ) )
  {
    return( 0 )
  }
  else
  {
    alleles <- table( unlist( lapply( seq( 1 , ncol( x ) ) , alleleSplit , as.character( x[i,] ) , ifelse1( is.null( Chr ) , NULL , Chr[i] ) , gender ) ) )
    maf <- min( alleles ) / sum( alleles )
    return( ifelse( is.nan( maf ) , NA , ifelse( isTRUE( maf < 1 ) , maf , 0 ) ) )
  }
}

#' Alternative ifelse
#'
#' Custom ifelse statement
#'
#' @param x boolean statement
#' @param a outcome if x is TRUE
#' @param b outcome if x is FALSE
#'
#' @return value of a or b
#'
#' @examples
#' \dontrun{
#' ifelse1(TRUE,1,0)==1
#' ifelse1(FALSE,1,0)==0
#' }
ifelse1 <- function( x , a , b )
{
  if ( x )
  {
    return( a )
  }
  else
  {
    return( b )
  }
}


#' Split Alleles for MAF
#'
#' Split Alleles
#'
#' @param i index of snp in X
#' @param x vector of genotypes for a single snp
#' @param Chr Chromosome each the snp is from
#' @param gender vector of the genders of each sample

#'
#' @return character vector of alleles
#'
#' @examples
#' \dontrun{
#'   exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=4,byrow=TRUE)
#'  exampleChr<-1
#'  exampleGender<-c("M","M","F","M","F","F","F","F","M")
#'
#'  identical(alleleSplit(1,exampleSNP[1,],exampleChr,exampleGender),c("A","A"))
#'  identical(alleleSplit(2,exampleSNP[1,],"X",exampleGender),c("A"))
#' }
alleleSplit <- function( i , x , Chr = NULL , gender = NULL )
{
  alleles <- strsplit( x[i] , "/" , fixed = TRUE )[[1]]
  if ( is.null( Chr ) )
  {
    return( alleles )
  }
  else
  {
    if ( Chr == "X" )
    {
      if ( gender[i] %in% c( "Male" , "male" , "m" , "M" ) )
      {
        return( alleles[1] )
      }
      else
      {
        return( alleles )
      }
    }
    else
    {
      return( alleles )
    }
  }
}

#' Calculate MAF for Illumina Data
#'
#' Calculate MAF for snps from Illumina Chips
#'
#' @param x vector of genotypes
#' @param Chr vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp

#'
#' @return numeric Minor allele frequency
#'
#' @examples
#' \dontrun{
#'  exampleSNP<-c(0,1,0,1,2,0,0,2,0,0,1,0,1,2,0,0,2,0,0,1,0,1,2,0,0,2,0)
#'
#'  mafFreqIllumina(exampleSNP,1,c("M","M","F","M","F","F","F","F","M"))
#' }
mafFreqIllumina <- function( x , Chr = NULL , gender = NULL )
{
  if ( isTRUE( all( is.na( x ) ) || all( x > 4 ) ) )
  {
    return( 0 )
  }
  else
  {
    x <- factor( x , levels = c( 0 , 1 , 2 ) , labels = c( "A/A" , "A/B" , "B/B" ) )
    alleles <- table( unlist( lapply( seq( 1 , length( x ) ) , alleleSplit , as.character( x ) , Chr , gender ) ) )
    maf <- min( alleles ) / sum( alleles )
    return( ifelse( is.nan( maf ) , NA , ifelse( isTRUE( maf < 1 ) , maf , 0 ) ) )
  }
}

#' Calculate MAF in parallel
#'
#' Helper function to fascilitate calculating MAF in parallel
#'
#' @param i index of snps in matrix to calculate MAF for
#' @param x matrix of genotypes
#' @param Chr vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp

#'
#' @return vector of named Minor Allele Frequencies
#'
#' @examples
#' \dontrun{
#'  exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,0,1,0,1,2,0,0,2,0,0,1,0,1,2,0,0,2,0),nrow=3,byrow=TRUE)
#'
#'  mafFreqIlluminaSF(c(1:3),exampleSNP,c(1,2,3),c("M","M","F","M","F","F","F","F","M"))
#' }
mafFreqIlluminaSF <- function( i , x , Chr = NULL , gender = NULL )
{
  return( unlist( apply( x[i,] , 1 , mafFreqIllumina , Chr , gender ) ) )
}

#' test Minor Allele Frequencies
#'
#' Calculate MAF for each snp, and identify SNPs failing to pass the Cutoff
#'
#' @param x Matrix of genotypes
#' @param qCutOff cutoff for Minor Allele Frequency
#' @param snowFall Boolean Value, should the MAF calculations be performed in parallel(TRUE) or series(FALSE)
#' @param Chr vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp
#'
#' @return vector of snps that failed to pass MAF test
#'
#' @examples
#' \dontrun{
#'   exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/A","A/A","A/A"),nrow=4,byrow=TRUE)
#'
#'  mafTest(exampleSNP, qCutOff = 0.5, snowFall=FALSE,Chr = c(1,1,2,2),gender = c("M","M","F","M","F","F","F","F","M","F"))
#' }
#' @export
mafTest <- function( x , qCutOff = 0.05 , snowFall = FALSE , Chr = NULL , gender = NULL )
{
  if ( inherits( x , c( "CNSet" , "SnpSuperSet" ) ) )
  {
    if ( snowFall )
    {
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "ff" , ( .packages() ) ) ) ) ) ) ) )
      {
        sfLibrary( ff , keep.source = FALSE )
      }
      fd <- calls( x )
      finalizer( fd )
      finalizer( fd ) <- "close"
      sfExport( "fd" )
      sfClusterEval( open( fd ) )
      sfExport( "alleleSplit" )
      sfExport( "mafFreqIllumina" )
      sfIndex <- unlist( sfLapply( chunk.ffdf( fd , by = nrow( fd ) / length( sfGetCluster() ) ) , mafFreqIlluminaSF , fd , Chr , gender ) )
      sfClusterEval( close( fd ) )

      return( which( sfIndex < qCutOff ) )
    }
    else
    {
      return( which( isTRUE( unlist( apply( calls( x )[] , 1 , mafFreqIllumina , Chr , gender ) ) < qCutOff ) ) )
    }
  }
  else
  {
    if ( snowFall )
    {
      return( which( unlist( sfLapply( seq( 1 , nrow( x ) ) , mafFreq , x , Chr , gender ) ) < qCutOff ) )
    }
    else
    {
      return( which( unlist( lapply( seq( 1 , nrow( x ) ) , mafFreq , x , Chr , gender ) ) < qCutOff ) )
    }
  }
}


#' Calculate Call rate
#'
#' Calculate call rates of the vector X, return true or false if meets call rate cutoff
#'
#' @param x vector of genotypes
#' @param crCutOff cutoff for callrate
#'
#' @return vector of snps that failed to pass Call rate test, True means failure
#'
#' @examples
#' \dontrun{
#'   examplePASS<-c("A/A","A/B","A/A","A/B",NA,"B/B","A/A","A/A","B/B","A/A","A/B","A/A")
#'   exampleFAIL<-c(NA,NA,NA,NA,NA,NA,"A/A","B/B","B/B","A/A","B/B","B/B")
#'
#'  callRate(examplePASS,.8)==FALSE
#'  callRate(exampleFAIL,.8)==TRUE
#' }
callRate <- function( x , crCutOff )
{
  return( ( 1 - sum( is.na( x ) ) / length( x ) ) <= crCutOff )
}

#' Calculate Call rate for Illumina
#'
#' Calculate call rates of the vector X, return true or false if meets call rate cutoff
#'
#' @param x vector of genotypes
#' @param crCutOff cutoff for callrate
#'
#' @return vector of snps that failed to pass Call rate test, True means failure
#'
#' @examples
#' \dontrun{
#'   examplePASS<-c(1,2,1,2,NA,3,1,1,3,1,2,1)
#'   exampleFAIL<-c(NA,NA,NA,NA,NA,NA,1,3,3,1,3,3)
#'
#'  callRateIllumina(examplePASS,.8)==FALSE
#'  callRateIllumina(exampleFAIL,.8)==TRUE
#' }
callRateIllumina <- function( x , crCutOff )
{
  return( ( 1 - sum( is.na( x ) ) / length( x ) ) <= crCutOff )
}

#' Calculate Call rate for Illumina -Helper func
#'
#' helper function to Calculate call rates of the vector X, return true or false if fails to meets call rate cutoff, via Snowfall
#'
#' @param i index in x to run call rate test on
#' @param x vector of genotypes
#' @param crCutOff cutoff for callrate
#'
#' @return vector of snps that failed to pass Call rate test, True means failure
#'
#' @examples
#' \dontrun{
#'   examplePASS<-c(1,2,1,2,NA,3,1,1,3,1,2,1)
#'   exampleFAIL<-c(NA,NA,NA,NA,NA,NA,1,3,3,1,3,3)
#'
#'  callRateIllumina(examplePASS,.8)==FALSE
#'  callRateIllumina(exampleFAIL,.8)==TRUE
#' }
callRateIlluminaSF <- function( i , x , crCutOff )
{
  return( unlist( apply( x[i,] , 1 , callRateIllumina , crCutOff ) ) )
}


#' Calculate Call rate -Master
#'
#' Master function to Calculate call rates of the snps in x, return true or false if fails to meets call rate cutoff
#'
#' @param x Matrix of genotypes
#' @param crCutOff cutoff for callrate
#' @param snowFall Boolean value. determine whether to calculate callrate in series(FALSE) or Parallel(TRUE)
#'
#' @return vector of snps that failed to pass Call rate test
#'
#' @examples
#' \dontrun{
#'   exampleSNP<-matrix(c("A/A","A/B","A/A",NA,"A/B","B/B","A/A","A/A","B/B",NA,NA,NA,
#'  "A/A","A/B","A/A",NA,"A/A","B/B","B/B","A/A",NA,"B/B","B/B","A/A","A/A",
#'  "B/B",NA,"A/A",NA,"A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",NA,NA,
#'  "A/A",NA,"A/B","B/B","A/A","A/A","B/B",NA,"A/A","A/A",NA,"A/A","A/A",NA,"A/A","A/A"),nrow=8,byrow=TRUE)
#'  rownames(exampleSNP)<-c("snp1","snp2","snp3","snp4","snp5","snp6","snp7","snp8")
#'
#'  all.equal(callRateTest(exampleSNP,.8),c("snp2"=2,"snp3"=3,"snp6"=6,"snp8"=8))
#' }
#'
#' @export
callRateTest <- function( x , crCutOff = 0.95 , snowFall = FALSE )
{
  if ( inherits( x , c( "CNSet" , "SnpSuperSet" ) ) )
  {
    if ( snowFall )
    {
      if ( any( packageLoaded <- unlist( sfClusterEval( !isTRUE( as.logical( grep( "ff" , ( .packages() ) ) ) ) ) ) ) )
      {
        sfLibrary( ff , keep.source = FALSE )
      }
      fd <- calls( x )
      finalizer( fd )
      finalizer( fd ) <- "close"
      sfExport( "fd" )
      sfClusterEval( open( fd ) )
      sfExport( "callRateIllumina" )
      sfIndex <- unlist( sfLapply( chunk.ffdf( fd , by = nrow( fd ) / length( sfGetCluster() ) ) , callRateIlluminaSF , fd , crCutOff ) )
      sfClusterEval( close( fd ) )
      return( which( sfIndex ) )
    }
    else
    {
      return( which( unlist( apply( calls( x )[] , 1 , callRateIllumina , crCutOff ) ) ) )
    }
  }
  else
  {
    if ( snowFall )
    {
      return( which( unlist( sfApply( x , 1 , callRate , crCutOff ) ) ) )
    }
    else
    {
      return( which( apply( x , 1 , callRate , crCutOff ) ) )
    }
  }
}
