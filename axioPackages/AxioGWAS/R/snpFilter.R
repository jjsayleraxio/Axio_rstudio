#' snp Filter
#'
#' Determine if snp passes Filtering Steps (MAF, CallRate), if so, return HWE pvalue
#'
#' @param i index of snp in x to apply snpFilter on
#' @param x Matrix of genotypes in "#/#" format
#' @param crCutOff Call Rate CutOff
#' @param mafCutOff Minor allele Frequency Cutoff
#' @param Chr vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp
#' 
#' @import ff
#'
#' @return dataframe of SNP index and HWE pvalue
#'
#' @examples
#' \dontrun{
#'   exampleGender<-c("Male","Male","Female","Female","Male","Male","Female","Female","Male","Male","Female")
#'   exampleSNP1<-c("A/A","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B")#callrate=100%,MAF=.4545, will return a HWE value
#'   exampleSNP2<-c("A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/B")#callrate=100%,MAF=.04545, will fail MAF
#'   exampleSNP3<-c(NA,NA,NA,"A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/B")#callrate=72.72%, will fail callRate#'   
#'   exampleSNPMatrix<-matrix(c(exampleSNP1,exampleSNP2,exampleSNP3),nrow=3,byrow=TRUE,dimnames=list(c("SNP1","SNP2","SNP3")))
#'   
#'   exampleChr<-c(1,1,1)
#'   
#'   snpFilter(1,exampleSNPMatrix,crCutOff=.95,mafCutOff=.05,Chr=exampleChr,gender=exampleGender)
#'   
#'   snpFilter(2,exampleSNPMatrix,crCutOff=.95,mafCutOff=.05,Chr=exampleChr,gender=exampleGender)
#'   
#'   snpFilter(3,exampleSNPMatrix,crCutOff=.95,mafCutOff=.05,Chr=exampleChr,gender=exampleGender)
#'   
#' }
#' @family snpFilter
#' @export
snpFilter<-function( i , x , crCutOff, mafCutOff , Chr , gender  ){
  UseMethod("snpFilter",x)
}

#' @family snpFilter
#' @export
snpFilter.default<-function( i , x , crCutOff, mafCutOff , Chr , gender ){
  stop(paste0("No Method for filtering snps for object of class ",class(x)[1],"."))
}

#' @family snpFilter
#' @export
snpFilter.matrix <- function( i , x , crCutOff, mafCutOff , Chr , gender  ){
  x<-x[i,]
  if(!callRate(x,crCutOff)){
    MAF<-mafFreq(1,matrix((x),nrow=1),Chr[i],gender)
    if(MAF>mafCutOff){
      return(data.frame( snp = i , MAF = MAF , HWE = hardyTest(x,FALSE) ))
    }
  }
}

#' @family snpFilter
#' @export
snpFilter.data.frame<-snpFilter.matrix
#' @family snpFilter
#' @export
snpFilter.data.table<-snpFilter.matrix

#' @family snpFilter
#' @export
snpFilter.ffdf <- snpFilter.dbDF <- function( i , x , crCutOff, mafCutOff , Chr , gender  )
{
  x<-as.character(unlist(c(x[i,])))
  if(!callRate(x,crCutOff)){
    MAF<-mafFreq(1,matrix(x,nrow=1),Chr[i],gender)
    if(MAF>mafCutOff){
      return(data.frame( snp = i , MAF = MAF , HWE = hardyTest(x,FALSE) ))
    }
  }
}

#' snp Filter Illumina
#'
#' Determine if snp passes Filtering Steps (MAF, CallRate), if so, return HWE pvalue, for Illumina
#'
#' @param i index of snp in x to apply snpFilter on
#' @param x Matrix of genotypes in "#/#" format
#' @param crCutOff Call Rate CutOff
#' @param mafCutOff Minor allele Frequency Cutoff
#' @param Chr vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp

#'
#' @return dataframe of SNP index and HWE pvalue
#'
#' @examples
#' \dontrun{
#'   exampleGender<-c("Male","Male","Female","Female","Male","Male","Female","Female","Male","Male","Female")
#'   exampleSNP1<-c(0,1,1,1,1,1,1,1,1,1,1)#callrate=100%,MAF=.4545, will return a HWE value
#'   exampleSNP2<-c(0,0,0,0,0,0,0,0,0,0,1)#callrate=100%,MAF=.04545, will fail MAF
#'   exampleSNP3<-c(NA,NA,NA,0,0,0,0,0,0,0,1)#callrate=72.72%, will fail callRate#'   
#'   exampleSNPMatrix<-matrix(c(exampleSNP1,exampleSNP2,exampleSNP3),nrow=3,byrow=TRUE,dimnames=list(c("SNP1","SNP2","SNP3")))
#'   
#'   exampleChr<-c(1,1,1)
#'   
#'   snpFilterIllumina(1,exampleSNPMatrix,crCutOff=.95,mafCutOff=.05,Chr=exampleChr,gender=exampleGender)
#'   
#'   snpFilterIllumina(2,exampleSNPMatrix,crCutOff=.95,mafCutOff=.05,Chr=exampleChr,gender=exampleGender)
#'   
#'   snpFilterIllumina(3,exampleSNPMatrix,crCutOff=.95,mafCutOff=.05,Chr=exampleChr,gender=exampleGender)
#'   
#' }
#' @export
snpFilterIllumina <- function( i , x , crCutOff, mafCutOff , Chr , gender  )
{
  x<-c(x[i,])
  if(!callRate(x,crCutOff)){
    MAF<-mafFreqIllumina(x,Chr[i],gender)
    if(MAF>mafCutOff){
      return(data.frame( snp = i , MAF = MAF , HWE = hardyTestIllumina(x,FALSE) ))
    }
  }
}

snpFilterIlluminaSF <- function( i , x , crCutOff, mafCutOff , Chr , gender  )
{
  return( do.call('rbind', lapply( i , snpFilterIllumina , x , crCutOff , mafCutOff , Chr , gender ) ) )
}

#' Run snpFilter
#'
#' Identify if a snp passes All Filtering Steps
#'
#' @param x vector of genotypes in "#/#" format
#' @param crCutOff Call Rate CutOff
#' @param mafCutOff Minor allele Frequency Cutoff
#' @param hweCutOff Hardy-Weinberg Equilibrium P-Value Cutoff
#' @param Chr named vector of which chromosome each single snp is in
#' @param gender vector of the genders of the sample for each snp
#' @param adjust pvalue adjustment method
#' @param snowFall Perform testing in parallel(TRUE) or in series(FALSE)
#' @param returnAll return HWE pvalues for all passing SNP
#'
#' @return vector of indeces of passing SNP
#'
#' @examples
#' \dontrun{
#'   exampleGender<-c("Male","Male","Female","Female","Male","Male","Female","Female","Male","Male","Female")
#'   
#'   exampleSNP1<-c("A/A","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B")#callrate=100%,MAF=.4545, will return a HWE value
#'   exampleSNP2<-c("A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/B")#callrate=100%,MAF=.04545, will fail MAF
#'   exampleSNP3<-c(NA,NA,NA,"A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/B")#callrate=72.72%, will fail callRate
#'   exampleSNP4<-c("A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B","A/B")#callrate=100%,MAF=.5, will return a HWE value
#'   exampleSNPMatrix<-matrix(c(exampleSNP1,exampleSNP2,exampleSNP3,exampleSNP4),nrow=4,byrow=TRUE,dimnames=list(c("SNP1","SNP2","SNP3","SNP4")))
#'   
#'   exampleChr<-c("SNP1"=1,"SNP2"=1,"SNP3"=1,"SNP4"=1)
#'   
#'   runsnpFilter(exampleSNPMatrix,crCutOff=.95,mafCutOff=.05, hweCutOff=0.05,Chr=exampleChr,gender=exampleGender,returnAll=TRUE)
#'   KeepSNPS<-runsnpFilter(exampleSNPMatrix,crCutOff=.95,mafCutOff=.05, hweCutOff=0.05,Chr=exampleChr,gender=exampleGender,returnAll=FALSE)
#'   
#'   #make exampleSNPMatrix into a SnpSuperSet
#'   exampleSNPMatrix<-matrix(as.numeric(gsub("A/B",1,gsub("A/A",0,exampleSNPMatrix))),nrow=4,byrow=FALSE,dimnames=list(c("SNP1","SNP2","SNP3","SNP4")))
#'   exampleSNPSS<-new("SnpSuperSet",alleleA=exampleSNPMatrix,alleleB=exampleSNPMatrix,call=exampleSNPMatrix,callProbability=exampleSNPMatrix)
#'   runsnpFilter(exampleSNPSS,crCutOff=.95,mafCutOff=.05, hweCutOff=0.05,Chr=exampleChr,gender=exampleGender,returnAll=TRUE)
#'   KeepSNPS<-runsnpFilter(exampleSNPSS,crCutOff=.95,mafCutOff=.05, hweCutOff=0.05,Chr=exampleChr,gender=exampleGender,returnAll=FALSE)
#'   
#'   
#' }
#' @family runsnpFilter
#' @export
runsnpFilter<-function(x, crCutOff, mafCutOff , hweCutOff , Chr  , gender, adjust = c( "holm" , "hochberg" , "hommel" , "bonferroni" , "BH" , "BY" , "fdr" , "none" ) , snowFall = FALSE , returnAll = FALSE){
  UseMethod("runsnpFilter",x)
}
#' @family runsnpFilter
#' @export
runsnpFilter.default<-function(x, crCutOff, mafCutOff , hweCutOff , Chr  , gender, adjust = c( "holm" , "hochberg" , "hommel" , "bonferroni" , "BH" , "BY" , "fdr" , "none" ) , snowFall = FALSE , returnAll = FALSE){
  if ( snowFall )
  {
    sfExport( "snpFilter" )
    snpFiltered<- do.call('rbind',sfLapply( rownames( x ) , snpFilter , x , crCutOff, mafCutOff , Chr , gender ) )
  }
  else
  {
    # return( p.adjust( unlist( apply( x , 1 , snpFilter , crCutOff, mafCutOff , Chr  , gender  ) ) , adjust ) )
    snpFiltered<- do.call('rbind',lapply( rownames( x )  , snpFilter , x , crCutOff, mafCutOff , Chr  , gender  ) )
  }
  
  snpFiltered$AdjHWE<-p.adjust(snpFiltered$HWE,adjust)
  if(!returnAll){
    return(which(as.character(rownames(x))%in%as.character(snpFiltered[which(snpFiltered$AdjHWE>hweCutOff),"snp"])))
  }else{
    return(snpFiltered)
  }
  
}
#' @family runsnpFilter
#' @export
runsnpFilter.CNSet <- function( x , crCutOff, mafCutOff , hweCutOff , Chr  , gender, adjust = c( "holm" , "hochberg" , "hommel" , "bonferroni" , "BH" , "BY" , "fdr" , "none" ) , snowFall = FALSE , returnAll = FALSE)
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
      sfExport( "snpFilterIllumina" )
      sfExport( "snpFilter" )
      snpFiltered <- do.call('rbind', sfLapply( bit::chunk( from = 1 , to = nrow( fd ) , by = nrow( fd ) / length( sfGetCluster() ) ) , snpFilterIlluminaSF , fd , crCutOff , mafCutOff, Chr , gender) )
      sfClusterEval( close( fd ) )
    }
    else
    {
      snpFiltered<- do.call('rbind',lapply( rownames( calls( x ) ) , snpFilterIllumina , calls( x ) , crCutOff, mafCutOff , Chr , gender ) )
    }
  
  snpFiltered$AdjHWE<-p.adjust(snpFiltered$HWE,adjust)
  if(!returnAll){
    return(which(as.character(rownames(x))%in%as.character(snpFiltered[which(snpFiltered$AdjHWE>hweCutOff),"snp"])))
  }else{
    return(snpFiltered)
  }
}
#' @family runsnpFilter
#' @export
runsnpFilter.SnpSuperSet <- runsnpFilter.CNSet


