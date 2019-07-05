#' Hardy-Weinberg Equilibrium
#'
#' Calculate the Hardy-Weinberg Equilibrium for a SNP
#'
#' @param x vector of genotypes for a single snp
#' @param returnAll Boolean value. If TRUE, returns the output of the fisher test. If FALSE, returns only the P-Values  from the fisher test.
#'
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=2,byrow=TRUE)
#'  
#'  hardyTest(exampleSNP[1,],TRUE)
#' }
#' @export
hardyTest <- function( x , returnAll = TRUE )
{
  if ( all( is.na( x ) ) )
  {
    return( 0 )
  }
  else
  {
    alleles <- table( unlist( lapply( as.character( x ) , function( x ) strsplit( x , "/" , fixed = TRUE ) ) ) )
    alleleSet <- FALSE
    if ( length( alleles ) < 2 )
    {
      if ( returnAll )
      {
        cat( "Monomorphic Allele" )
        print( alleles )
        return( 1 )
      }
      else
      {
        return( 1 )
      }
    }
    else
    {
      genotypes <- table( as.character(x) )
      if ( isTRUE( all.equal( as.integer( alleles[1] ) , as.integer( alleles[2] ) ) ) )
      {
        
        major <- alleles[1]
        minor <- alleles[2]
      }
      else
      {
        major <- alleles[which.max( alleles )]
        minor <- alleles[which.min( alleles )]
      }
      p <- major / sum( alleles )
      q <- 1 - p
      het <- c( paste( names( major ) , names( minor ) , sep = "/" ) , paste( names( minor ) , names( major ) , sep = "/" ) )
      het <- het[which( het %in% names( genotypes ) )]
      homwt <- genotypes[paste( names( major ) , names( major ) , sep = "/" )]
      het <- genotypes[het]
      if ( length( het ) < 1 )
      {
        het <- NA
      }
      hommt <- genotypes[paste( names( minor ) , names( minor ) , sep = "/" )]
      observed <- c( ifelse( is.na( homwt ) , 0 , homwt ) , ifelse( is.na( het ) , 0 , het ) , ifelse( is.na( hommt ) , 0 , hommt ) )
      expected <- round( c( p^2 , 2 * p * q , q^2 ) * sum( genotypes ) )
      freqs <- rbind( observed , expected )
      rownames( freqs ) <- c(  "Observed" , "Expected" )
      colnames( freqs ) <- c( paste( names( alleles )[which.max( alleles )] , names( alleles )[which.max( alleles )] , sep = "/" ) , paste( names( alleles )[which.max( alleles )] , names( alleles )[which.min( alleles )] , sep = "/" ) , paste( names( alleles )[which.min( alleles )] , names( alleles )[which.min( alleles )] , sep = "/" ) )
      if ( sum( freqs ) < 1 )
      {
        return( 1 )
      }
      if ( returnAll )
      {
        print( freqs )
        return( fisher.test( freqs ) )
      }
      else
      {
        return( fisher.test( freqs )$p.value )
      }
    }
  }
}

#' Hardy-Weinberg Equilibrium for Illumina Calls
#'
#' Calculate the Hardy-Weinberg Equilibrium for a SNP from an Illumina Chip
#'
#' @param x Vector of numeric (0,1,2) genotypes for a single SNP, 
#' @param returnAll Boolean value. If TRUE, returns the output of the fisher test. If FALSE, returns only the P-Values  from the fisher test.
#'
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,1,0,0,2,2,0,2,2,0,0,2,0,1,0,0,2,2,0,2,0,1,0,1,2,0,0,2),nrow=2,byrow=TRUE)
#'  
#'  hardyTestIllumina(exampleSNP[1,],TRUE)
#' }
hardyTestIllumina <- function( x , returnAll = TRUE )
{
  if ( all( is.na( x ) ) )
  {
    return( 0 )
  }
  else
  {
    x <- factor( x , levels = c( 0 , 1 , 2 ) , labels = c( "A/A" , "A/B" , "B/B" ) )
    return( hardyTest( x , returnAll ) )
  }
}

#' Hardy-Weinberg Equilibrium for Illumina Calls
#'
#' Calculate the Hardy-Weinberg Equilibrium for a SNP from an Illumina Chip
#'
#' @param x Vector of numeric (1,2,3) genotypes for a single SNP, 
#' @param returnAll Boolean value. If TRUE, returns the output of the fisher test. If FALSE, returns only the P-Values  from the fisher test.
#'
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,1,0,0,2,2,0,2,2,0,0,2,0,1,0,0,2,2,0,2,0,1,0,1,2,0,0,2),nrow=2,byrow=TRUE)
#'  
#'  hardyTestIllumina(exampleSNP[1,],TRUE)
#' }
hardyTestIlluminaSF <- function( i , x , returnAll = FALSE )
{
  return( unlist( apply( x[i,] , 1 , hardyTestIllumina , returnAll ) ) )
}

#' Hardy-Weinberg Equilibrium for ALL snps
#'
#' Calculate the Hardy-Weinberg Equilibrium for each SNP, and return adjusted P-Values
#'
#' @param x Matrix of genotypes for a
#' @param adjust The method to adjust the P-values by after finding the HWE for each SNP. See \code{\link[stats]{p.adjust}} for more options
#' @param snowFall Boolean Value, should the calculations be performed in parallel(TRUE) or series(FALSE)
#'
#' @return vector of p-values
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(
#'  c("A/A","A/A","A/A","A/A","A/A","A/A","A/A",
#'  "A/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A",
#'  "B/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A",
#'  "A/A","B/B","B/B","A/A","B/B","A/A","A/B","A/A",
#'  "A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A",
#'  "A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A",
#'  "A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A",
#'  "B/B","A/A","A/B","B/B","A/A","A/A","A/B","B/B","A/A",
#'  ),nrow=8,byrow=TRUE)
#'  
#'  hweTest(exampleSNP)
#' }
#' @export
hweTest <- function( x , adjust = c( "holm" , "hochberg" , "hommel" , "bonferroni" , "BH" , "BY" , "fdr" , "none" ) , snowFall = FALSE )
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
      sfExport( "hardyTestIllumina" )
      sfExport( "hardyTest" )
      sfPVal <- unlist( sfLapply( bit::chunk( from = 1 , to = nrow( fd ) , by = nrow( fd ) / length( sfGetCluster() ) ) , hardyTestIlluminaSF , fd ) )
      sfClusterEval( close( fd ) )
      return( p.adjust( sfPVal , adjust ) )
    }
    else
    {
      return( p.adjust( unlist( apply( calls( x )[] , 1 , hardyTestIllumina , FALSE  ) ) , adjust ) )
    }
  }
  else
  {
    if ( snowFall )
    {
      return( p.adjust( unlist( sfApply( x , 1 , hardyTest , FALSE ) ) , adjust ) )
    }
    else
    {
      return( p.adjust( unlist( apply( x , 1 , hardyTest , FALSE ) ) , adjust ) )
    }
  }
}

#' Monomorphic SNP
#'
#' Calculate if the SNP is monomorphic or not
#'
#' @param x vector of genotypes for a single snp
#' 
#' @return boolean value
#' 
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=2,byrow=TRUE)
#'  
#'  monomorphicSNP(exampleSNP[1,])
#' }
monomorphicSNP <- function( x )
{
  return( length( table( unlist( lapply( as.character( x ) , function( x ) strsplit( x , "/" , fixed = TRUE ) ) ) ) ) < 2 )
}

#' Monomorphic SNP for Illumina
#'
#' Calculate if the SNP is monomorphic or not from Illumina Chips
#'
#' @param x vector of genotypes for a single snp
#' 
#' @return boolean value
#' 
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,1,0,0,2,2,0,2,2,0,0,2,0,1,0,0,2,2,0,2,0,1,0,1,2,0,0,2),nrow=2,byrow=TRUE)
#'  
#'  monomorphicSNPIllumina(exampleSNP[1,])
#' }
monomorphicSNPIllumina <- function( x )
{
  return( length( table( x ) ) < 2 )
}

#' Monomorphic SNP for Illumina, Snowfall
#'
#' Helper function to calculate if the SNP is monomorphic or not from Illumina Chips, allow chunking by SnowFall package
#'
#' @param i index of Matrix to use in apply
#' @param x Matrix of genotypes, each row is a SNP, column is a sample
#' 
#' @return vector of named boolean values
#' 
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,1,0,0,2,2,0,2,2,0,0,2,0,1,0,0,2,2,0,2,0,1,0,1,2,0,0,2),nrow=6,byrow=TRUE)
#'  rownames(exampleSNP)<-c("snp1","snp2","snp3","snp4","snp5","snp6")
#'  
#'  monomorphicSNPIlluminaSF(c(1:6),exampleSNP)
#' }
monomorphicSNPIlluminaSF <- function( i , x )
{
  return( unlist( apply( x[i,] , 1 , monomorphicSNPIllumina ) ) )
}

#' Monomorphic SNP for Each SNP
#'
#' Function to calculate if the SNP is monomorphic or not for each SNP in x
#'
#' @param x Matrix of genotypes
#' @param snowFall Boolean, should monomorphic calculations be performed in parallel(TRUE) or series(FALSE).
#' 
#' @return named vector of monomorphic snps with index in x
#' 
#' @examples 
#' \dontrun{
#'  exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,1,0,0,2,2,0,2,2,0,0,2,0,1,0,0,2,2,0,2,0,1,0,1,2,0,0,2,0,0,0,0,0,0),nrow=7,byrow=TRUE)
#'  rownames(exampleSNP)<-c("snp1","snp2","snp3","snp4","snp5","snp6","snp7")
#'  
#'  findMonomorphic(exampleSNP)
#' }
findMonomorphic <- function( x , snowFall = FALSE )
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
      sfExport( "monomorphicSNPIllumina" )
      sfIndex <- unlist( sfLapply( chunk.ffdf( fd , BATCHBYTES = 2^14 ) , monomorphicSNPIlluminaSF , fd ) )
      sfClusterEval( close( fd ) )
      return( which( sfIndex ) )
    }
    else
    {
      return( which( unlist( apply( calls( x ) , 1 , monomorphicSNPIllumina ) ) ) )
    }
  }
  else
  {
    if ( snowFall )
    {
      return( which( unlist( sfApply( x , 1 , monomorphicSNP ) ) ) )
    }
    else
    {
      return( which( unlist( apply( x , 1 , monomorphicSNP ) ) ) )
    }
  }
}
