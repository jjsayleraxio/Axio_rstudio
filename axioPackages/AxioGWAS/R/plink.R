
#' Create a PED file
#'
#' create a PED file for Plink
#'
#' @param file the filename of the PED file
#' @param phenotype the categorical phenotype that is a column name in demographics
#' @param demographics data.frame containing phenotypes of the subjects. Samples by row, Phenotypes by column
#' @param genotypes an Object containing the genotypes. If not a CNset or SnpSuperSet, a matrix of Snp by row, subject by column
#'
#' @examples 
#' \dontrun{
#'  exampleDemographic<-data.frame(Phenotype1=c("Level1","Level1","Level2","Level2","Level1","Level2","Level1","Level1","Level2"),
#'                                SEX=c("F","M","M","F","F","F","M","M","M"),
#'                                row.names=c("Sub1","Sub2","Sub3","Sub4","Sub5","Sub6","Sub7","Sub8","Sub9"))
#'                                
#'  exampleGenotype<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B",
#'  "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","B/B","A/A","A/A",
#'  "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B",
#'  "A/A","A/B","B/B","A/A","A/A","B/B"),nrow=9,byrow=TRUE)
#'  rownames(exampleGenotype)<-c("Sub1","Sub2","Sub3","Sub4","Sub5","Sub6","Sub7","Sub8","Sub9")
#'  
#'  makePedFile(file="test",phenotype="Phenotype1",demographics=exampleDemographic,genotypes=exampleGenotype)
#'  
#' }
makePedFile <- function( file = "ped" , phenotype = "ACSTRESN" , demographics , genotypes )
{
  file <- paste( paste( file , phenotype , sep = "-" ) , "ped" , sep = "." )
  rNames <- rownames( demographics )
  if ( inherits( genotypes , c( "CNSet" , "SnpSuperSet" ) ) )
  {
    for ( i in seq( 1 , nrow( demographics ) ) )
    {
      gen <- unlist( lapply( calls( genotypes[,match( rNames[i] , sampleNames( genotypes ) )] ) , plinkifyGenotypes ) )
      cat( c( i , rNames[i] , "0" , "0" , ifelse( demographics[i,"SEX"] == "M" , 1 , 2 ) , demographics[i,phenotype] , gen , "\n" ) , file = file , sep = "\t" , append = TRUE )
    }
  }
  else
  {
    genotypes <- genotypes[rNames,]
    
    ped <- as.data.frame( matrix( nrow = length( rNames ) , ncol = 6 + ncol( genotypes ) ) )
    for ( i in seq( 1 , nrow( demographics ) ) )
    {
      gen <- gsub( "/" , " " , unlist( lapply( genotypes[i,] , as.character ) ) )
      ped[i,] <- c( i , rNames[i] , "0" , "0" , ifelse( demographics[i,"SEX"] == "M" , 1 , 2 ) , demographics[i,phenotype] , gen )
    }
    write.table( file = file , ped , sep = "\t" , row.names = FALSE , quote = FALSE , col.names = FALSE , na = "0 0" )
  }
}


#' Plink-ify genes
#'
#' convert numeric genotypes to plink format
#'
#' @param x numeric genotype
#'
#' @examples 
#' \dontrun{
#'  exampleSNP<-c(1,2,1,2,3,1,1,3,1)
#'  
#'  plinkifyGenotypes(exampleSNP[1])=="A A"  
#' }
plinkifyGenotypes <- function( x )
{
  if ( is.na( x ) )
  {
    return( "0 0" )
  }
  else
  {
    if ( x == 1 )
    {
      return( "A A" )
    }
    else if ( x == 2 )
    {
      return( "A B" )
    }
    else
    {
      return( "B B" )
    }
  }
}

#' Create an Alternate Pheno file
#'
#' create a Phenotype file for Plink
#'
#' @param file the filename of the Phenotype file (will append a .phe to the end)
#' @param phenotypes vector of phenotype names that is a column name in demographics to keep
#' @param demographics data.frame containing phenotypes of the subjects. Samples by row, Phenotypes by column
#'
#' @examples 
#' \dontrun{
#'  exampleDemographic<-data.frame(Phenotype1=c("Level1","Level1","Level2","Level2","Level1","Level2","Level1","Level1","Level2"),
#'                                 Phenotype2=c(1,1,2,2,5,-1,2,-6,12),
#'                                 Phenotype3=c("Level5","Leve2","Level2","Level3","Level10","Level1","Level7","Level1","Level1"),
#'                                SEX=c("F","M","M","F","F","F","M","M","M"),
#'                                row.names=c("Sub1","Sub2","Sub3","Sub4","Sub5","Sub6","Sub7","Sub8","Sub9"))
#'  
#'  makeAltPhenFile(file="test",phenotypes=c("Phenotype1","Phenotype2","Phenotype3"),demographics=exampleDemographic)
#'  
#' }
makeAltPhenFile <- function( file = "ped" , phenotypes = c( "RESPSTD" , "RESPBAD" ) , demographics )
{
  file <- paste( file , "phe" , sep = "." )
  rNames <- rownames( demographics )
  ped <- as.data.frame( matrix( nrow = length( rNames ) , ncol = 2 + length( phenotypes ) ) )
  for ( i in phenotypes )
  {
    if ( inherits( demographics[,i] , "factor" ) )
    {
      demographics[,i] <- unlist( gsub( " " , "_" , as.character( demographics[,i] ) , fixed = TRUE ) )
    }
  }
  for ( i in seq( 1 , nrow( demographics ) ) )
  {
    ped[i,] <- c( i , rNames[i] , demographics[i,phenotypes] )
  }
  names( ped ) <- c( "FID" , "IID" , phenotypes )
  write.table( file = file , ped , sep = "\t" , row.names = FALSE , quote = FALSE , na = "-9" )
}


#' Create an Covariate file
#'
#' create the covariate file for Plink
#'
#' @param file the filename of the Covariates file (will append a .cvr to the end)
#' @param covariates vector of phenotype names that is a column name in demographics to add to the covariates file
#' @param demographics data.frame containing phenotypes of the subjects. Samples by row, Phenotypes by column
#'
#' @examples 
#' \dontrun{
#'  exampleDemographic<-data.frame(Phenotype1=c("Level1","Level1","Level2","Level2","Level1","Level2","Level1","Level1","Level2"),
#'                                 Phenotype2=c(1,1,2,2,5,-1,2,-6,12),
#'                                 Phenotype3=c("Level5","Leve2","Level2","Level3","Level10","Level1","Level7","Level1","Level1"),
#'                                SEX=c("F","M","M","F","F","F","M","M","M"),
#'                                row.names=c("Sub1","Sub2","Sub3","Sub4","Sub5","Sub6","Sub7","Sub8","Sub9"))
#'  
#'  makeAltPhenFile(file="test",covariates=c("Phenotype1","Phenotype2","Phenotype3"),demographics=exampleDemographic)
#'  
#' }
makeCovarFile <- function( file = "ped" , covariates = c( "GENDER" ) , demographics )
{
  file <- paste( file , "cvr" , sep = "." )
  rNames <- rownames( demographics )
  ped <- as.data.frame( matrix( nrow = length( rNames ) , ncol = 2 + length( covariates ) ) )
  for ( i in covariates )
  {
    if ( inherits( demographics[,covariates] , c( "factor" , "character" ) ) )
    {
      demographics[,i] <- as.character( unlist( gsub( " " , "_" , as.character( demographics[,i] ) , fixed = TRUE ) ) )
    }
  }
  for ( i in seq( 1 , nrow( demographics ) ) )
  {
    ped[i,] <- c( i , rNames[i] , unlist( lapply( demographics[i,covariates] , dumpFactorToCharacter ) ) )
  }
  names( ped ) <- c( "FID" , "IID" , covariates )
  write.table( file = file , ped , sep = "\t" , row.names = FALSE , quote = FALSE , na = "-9" )
}
#' Convert Factor to character
#'
#' convert factors to character for the covar file
#'
#' @param x A single character or factor instance
#'
#' @examples 
#' \dontrun{
#'  exampleDemographic<-data.frame(Phenotype1=c("Level1","Level1","Level2","Level2","Level1","Level2","Level1","Level1","Level2"),
#'                                 Phenotype2=c(1,1,2,2,5,-1,2,-6,12),
#'                                 Phenotype3=c("Level5","Leve2","Level2","Level3","Level10","Level1","Level7","Level1","Level1"),
#'                                SEX=c("F","M","M","F","F","F","M","M","M"),
#'                                row.names=c("Sub1","Sub2","Sub3","Sub4","Sub5","Sub6","Sub7","Sub8","Sub9"))
#'  
#'  identical(dumpFactorToCharacter(exampleDemographic["Sub1","Phenotype3"]),"Level5")
#' }
dumpFactorToCharacter <- function( x )
{
  return( ifelse( inherits( x , c( "factor" , "character" ) ) , unlist( gsub( " " , "_" , as.character( x ) , fixed = TRUE ) ) , x ) )
}


#' Create MAP file for Plink
#'
#' Converts SNP data and phenotype to MAP file for use in Plink
#'
#' @param file Name of the MAP file, will end up being 'file'-phenotype.map
#' @param phenotype The name of the phenotype to include in the map file name
#' @param snpInfo A data.frame of data about the SNPs, including a chr column and a pos column. snp names are rownames
#'
#' @examples 
#' \dontrun{
#'  exampleSNPInfo<-data.frame(chr=c("1","1","2","2","3","3","4","4","5","5"),
#'                             pos=c(1,2,1,2,1,2,1,2,1,2),
#'                             row.names=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10"))
#'  
#'  makeMapFile(file="test",phenotype="Phenotype1",snpInfo=exampleSNPInfo)
#' }
makeMapFile <- function( file = "ped" , phenotype = "ACSTRESN" , snpInfo )
{
  file <- paste( paste( file , phenotype , sep = "-" ) , "map" , sep = "." )
  if ( inherits( snpInfo , c( "CNSet" , "SnpSuperSet" ) ) )
  {
    snpinfo <- data.frame( chromosome( snpInfo ) , featureNames( snpInfo ) , position( snpInfo ) )
  }else{
    snpinfo <- data.frame( gsub( "chr" , "" , snpInfo$chr ) , rownames( snpInfo ) , snpInfo$pos )
  }
  write.table( file = file , snpinfo , sep = "\t" , row.names = FALSE , quote = FALSE , col.names = FALSE )
}
