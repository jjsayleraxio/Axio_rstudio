
#' Plot Results of 
#'
#' Converts SNP data and phenotype to MAP file for use in Plink
#'
#' @param x output from runGWAS of Pvalues
#' @param phenotype The name of the phenotype to include in the map file name
#' @param snpInfo A data.frame of data about the SNPs, including a chr column and a pos column. snp names are rownames
#'
#' @examples 
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,2,2,1,2,2),var2=c(2,2,2,3,3,3,2,2,3,2,3),var3=c(4,3,3,4,4,4,3,2,4,2,3),gender=c("Male","Female","Male","Female","Male","Female","Male","Male","Female","Female","Female"))
#' exampleSNP<-matrix(c("A/B","A/A","A/A","A/B","B/B","B/B","A/A","A/B","A/B","A/A","A/A",
#'                      "A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B","A/A",
#'                      "A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/B","A/B","B/B","A/A",
#'                      "B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B"),
#'                      nrow=4,byrow=TRUE,dimnames=list(c("SNP1","SNP2","SNP3","SNP4"),NULL))
#' exampleSNPInfo<-data.frame(Chr=c(1,2,3,4))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' exampleX<-runGWAS(pgx=examplePGx,snp=exampleSNP,Model=exampleModel,func = AxioGWAS::gwasFunc,interactionTerm="var2",genderVar="gender",databasePath="runGWAS.DB", dbTable="runGWAS_Test", snpInfo=exampleSNPInfo, overwriteEntry=TRUE, adjust="none")
#' 
#' plotResults(x=exampleX,pgx=examplePGx,genotypes=exampleSNP,interactionTerm="var2", snpInfo = exampleSNPInfo, genderVar = "gender")
#' }
plotResults <- function( x , pgx , genotypes , interactionTerm , pCutOff = 0.05 , snpInfo, genderVar)
{
  sigIndx <- which( x$adjPVal$adjpval < pCutOff )
  # numericGenotype<-lapply(sigIndx,function())
  summaryDF<-data.frame(do.call('rbind',lapply(sigIndx,function(indx,genotype,chromosome,gender,snpName,pgx,x){
    ADDVec<-data.frame(Add=c(makeAdditive( genotype[indx,] , chromosome[indx] ,gender)))
    sigDF <- data.frame( pgx , ADDVec )
    respVar <- rownames( attr( terms( formula( x$modelFits[[indx]]$lm ) ) , "factors" ) )[attr( terms( formula( x$modelFits[[indx]]$lm ) ) , "response" )]
    sigDF$Mean <- unlist( lapply( split( sigDF[,respVar] , interaction( sigDF[,interactionTerm] , sigDF$Add ) ) , function(x){return(mean(x,na.rm=TRUE))} ) )[paste( sigDF[,interactionTerm] , sigDF$Add , sep = "." )]
    sigDF$CI.half <- NA
    sigDF[names( predict( x$modelFits[[indx]]$lm , se.fit = TRUE )$fit ),"CI.half"] <- predict( x$modelFits[[indx]]$lm , se.fit = TRUE )$se / 2
    summaryDF <- unique( sigDF[!is.na( sigDF$Mean ),c( interactionTerm , "Add" , "Mean" , "CI.half" )] )
    names( summaryDF ) <- c( interactionTerm , "Additive_Effect" , "Mean" , "Half CI" )
    summaryDF[,"Additive_Effect"] <- as.factor( summaryDF[,"Additive_Effect"] )
    
    summaryDF$SNPName<-snpName[indx]
    return(summaryDF)
    },genotypes,snpInfo$Chr,as.character(pgx[,genderVar]),rownames(x$adjPVals),pgx,x)))
  
  eval(parse(text=paste0("plotRes<-qplot( data = summaryDF , y = Mean  , x = ",interactionTerm," , group = Additive_Effect , geom = c( \"point\" , \"line\" ) , color = Additive_Effect , shape = Additive_Effect , width = 0.25 )")))
  return(plotRes+facet_grid(.~SNPName))
}




plotCallConfidenceLog2 <- function( gtSet , snpID )
{
  invisible( open( gtSet ) )
  i <- match( snpID , featureNames( gtSet ) )
  genotypeSet <- gtSet[i,]
  invisible( close( gtSet ) )
  a <- as.matrix( log2( A( genotypeSet ) ) )
  b <- as.matrix( log2( B( genotypeSet ) ) )
  gt <- as.integer( calls( genotypeSet ) )
  gt.conf <- as.numeric( confs( genotypeSet ) )
  min.conf <- min( gt.conf )
  max.conf <- max( gt.conf )
  sc <- ( gt.conf - min.conf ) / ( max.conf - min.conf )
  gtData <- data.frame( a = t( a ) , b = t( b ) , gt = gt , sc = sc )
  names( gtData ) <- c( "a" , "b" , "gt" , "sc" )
  gtData$gt <- factor( gtData$gt )
  return( qplot( x = a , y = b , data = gtData , colour = gt , alpha = sc , main = substitute( "Plot of Genotype Call and Confidence Score by Alelle and" ~ Log[2] ~ "Intensity for SNP" ~ snpID , list( snpID = snpID) ) , xlab = expression( Log[2]~"Allele A Intensity" ) , ylab = expression( Log[2]~"Allele B Intensity" ) ) + labs( alpha = "Confidence\nScore" , colour = "Genotype" ) )
}

plotSNR <- function( x )
{
  invisible( open( x$SNR ) )
  snr <- x$SNR[]
  invisible( close( x$SNR ) )
  snr <- data.frame( snr = snr , Class = factor( ifelse( snr < 25 , "Bad" , "Good" ) ) )
  return( qplot( snr , data = snr , fill = Class , xlab = "Signal to Noise Ratio" , ylab = "Count" , main = "Histogram of Signal to Noise Ratio,\n< 25 is Indicator of Poor Sample Quality" , geom = c( "histogram" ) ) )
}


#' Plot Results of 
#'
#' Converts SNP data and phenotype to MAP file for use in Plink
#'
#' @param x dataframe containing chromosome and position for the snps intended to use
#' @param scale The name of the phenotype to include in the map file name
#' @param reference Which genome to use as reference
#' @param excludeXYM Boolean value, whether to keep (FALSE) or remove (TRUE) the X,Y and M Chromosomes
#'
#' @import hgu95av2.db
#'
#' @examples 
#' \dontrun{
#' exampleX<-data.frame(CHR=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,5),
#'                      pos=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,10,11,12,1))
#' unscaledXAxis<-scaleXAxis(x=exampleX,scale=FALSE,reference="hgu95av2",excludeXYM=TRUE)
#' scaledXAxis<-scaleXAxis(x=exampleX,scale=TRUE,reference="hgu95av2",excludeXYM=TRUE)
#' }
scaleXAxis <- function( x , scale = FALSE , reference = "hgu95av2" , excludeXYM = TRUE )
{
  if ( scale )
  {
    cLengths <- c( 0 , cumsum( table( x$CHR ) ) )
    cLengthNames <- names( cLengths )
    names( cLengths ) <- c( cLengthNames[-1] , "End" )
  }
  else
  {
    cLengths <- as.numeric( eval( parse( text = paste( reference , "CHRLENGTHS" , sep = "" ) ) ) )
    names( cLengths ) <- names( eval( parse( text = paste( reference , "CHRLENGTHS" , sep = "" ) ) ) )
    cLengths <- cLengths[-grep( "_" , names( cLengths ) , fixed = TRUE )]
    cLengths <- cLengths[order( factor( names( cLengths ) , levels = c( 1:22 , "X" , "Y" , "M" ) , ordered = TRUE ) )]
    cNames <- names( cLengths )
    cLengths <- c( 0 , cumsum( cLengths ) )
    names( cLengths ) <- c( cNames , "End" )
  }
  if ( excludeXYM )
  {
    if ( any( names( cLengths ) %in% c( "X" , "Y" , "M" ) ) )
    {
      cLengths <- cLengths[-( which( names( cLengths ) %in% c( "X" , "Y" , "M" ) ) + 1 )]
      names( cLengths )[length( cLengths )] <- "End"
    }
  }
  xLabelLocs <- vector( length = length( cLengths ) - 1 )
  for ( i in seq( 1 , length( xLabelLocs ) ) )
  {
    xLabelLocs[i] <- median( c( cLengths[i] , cLengths[i + 1] ) )
  }
  names( xLabelLocs ) <- names( cLengths[-length( cLengths )] )
  cLengths <- cLengths[sort( unique( x$CHR ) )]
  xLabelLocs <<- xLabelLocs[seq( 1 , max( which( c( 1:22 , "X" , "Y" , "M" ) %in% names( cLengths ) ) ) )]
  return( list( cLengths = cLengths , xLabelLocs = xLabelLocs ))
}


#' Plot Results of 
#'
#' Converts SNP data and phenotype to MAP file for use in Plink
#'
#' @param x a dataframe containing the columns: CHR, pos, pValue. the rownames are the snpname
#' @param type ?
#' @param chr vector of chromosomes to keep
#'
#' @examples 
#' \dontrun{

#' exampleX<-data.frame(CHR=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,5),
#'                      pos=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,10,11,12,1))
#' exampleX$pValue<-runif(dim(exampleX)[1])
#' manhattanPlot(x=exampleX,chr=c(1,2,3,4))
#' }
manhattanPlot <- function( x , type = "all" , chr = NULL , xlab = "Chromosome" , ylab = expression( -Log[10]~"p-Value" ) , main = "Manhattan Plot" , reference = "hgu95av2" , excludeXYM = TRUE , pCutOff = 0.05 , sCutOff = 0.01 , regions = NULL , rawChromosomes = FALSE , bleach = 0 )
{
  x$Start <- as.numeric( as.character( x[,"pos"] ) )
  x$pValue <- -log10( x$pValue )
  pCutOff <- -log10( pCutOff )
  xAxisInfo <- scaleXAxis( x , reference = reference , excludeXYM = excludeXYM )
  if ( length( ifelse( is.null( chr ) , runif( 2 ) , chr ) ) > 0 )
  {
    if ( excludeXYM )
    {
      for ( i in sort( unique( factor( x$CHR , levels = c( 1:22 ) , ordered = TRUE ) ) ) )
      {
        x$Start[which( x$CHR == i )] <- as.numeric( as.character( x$pos[which( x$CHR == i )] ) ) + xAxisInfo$cLengths[i]
      }
    }
    else
    {
      for ( i in sort( unique( factor( x$CHR , levels = c( 1:22 , "X" , "Y" , "M" ) , ordered = TRUE ) ) ) )
      {
        x$Start[which( x$CHR == i )] <- as.numeric( as.character( x$pos[which( x$CHR == i )] ) ) + xAxisInfo$cLengths[i]
      }
    }
  }
  x$Group <- NA
  x$Group[which( as.numeric( x$CHR ) %% 2 == 0 )] <- "red"
  x$Group[which( as.numeric( x$CHR ) %% 2 != 0 )] <- "black"
  if ( !is.null( chr ) )
  {
    chrMin <- min( xAxisInfo$cLengths[which( names( xAxisInfo$cLengths ) %in% chr )] )
    chrMax <- ifelse( names( xAxisInfo$cLengths )[max( which( names( xAxisInfo$cLengths ) %in% chr ) )] == names( xAxisInfo$cLengths )[length( xAxisInfo$cLengths )] , 2 * ( xAxisInfo$xLabelLocs[names( xAxisInfo$cLengths )[max( which( names( xAxisInfo$cLengths ) %in% chr ) )]] - chrMin ) , xAxisInfo$cLengths[max( which( names( xAxisInfo$cLengths ) %in% chr ) ) + 1] - chrMin )
    x <- x[which( x$CHR %in% chr ),]
    xAxisInfo$cLengths <- xAxisInfo$cLengths[which( names( xAxisInfo$cLengths ) %in% chr )] - chrMin
    xAxisInfo$xLabelLocs <- xAxisInfo$xLabelLocs[which( names( xAxisInfo$xLabelLocs ) %in% chr )] - chrMin
    regions <- regions[which( regions[,"chromosome_name"] %in% chr ),]
    x$Start <- x$Start - chrMin
  }
  sIndex <- x$pValue > pCutOff
  nsIndex <- x$pValue <= pCutOff
  if ( nrow( x ) < 10 )
  {
    hexPVals <- qplot( Start , pValue , data = x[which( nsIndex ),] , geom = "point" , position = "jitter" , main = main , xlab = xlab , ylab = ylab ) + scale_fill_continuous( "N" , low = "lightblue" , high = "darkblue" )
  }
  else
  {
    hexPVals <- qplot( Start , pValue , data = x[which( nsIndex ),] , geom = "hex" , main = main , xlab = xlab , ylab = ylab ) + scale_fill_continuous( "N" , low = "lightblue" , high = "darkblue" )
  }
  if ( length( which( sIndex ) ) > 0 )
  {
    hexPVals <- hexPVals + geom_point( aes( x = Start , y = pValue , colour = Group ) , data = x[which( sIndex ),] )
  }
  if ( !is.null( chr ) )
  {
    hexPVals <- hexPVals + scale_x_continuous( breaks = as.numeric( xAxisInfo$xLabelLocs ) , labels = as.character( names( xAxisInfo$xLabelLocs ) ) , limits = c( 0 , chrMax ) ) + scale_colour_discrete( guide = FALSE )
  }
  else
  {
    hexPVals <- hexPVals + scale_x_continuous( breaks = as.numeric( xAxisInfo$xLabelLocs ) , labels = as.character( names( xAxisInfo$xLabelLocs ) ) ) + scale_colour_discrete( guide = FALSE )
  }
  hexPVals <- hexPVals + geom_hline( yintercept = pCutOff , colour = "yellow" ) + geom_hline( yintercept = -log10( sCutOff ) , colour = "red" )
  hexPVals <- hexPVals + geom_text( aes( label = paste( "p < " , 10^-pCutOff ) , x = 1000 , y = pCutOff ) , size = 4 , vjust = 0 ) + geom_text( aes( label = paste( "p < " , sCutOff ) , x = 1000 , y = -log10( sCutOff ) ) , size = 4 , vjust = 0 )
  if ( !is.null( regions ) )
  {
    ymax <- 0.95 * max( x$pValue )
    regions$xmin <- regions$xmax <- regions$xmean <-regions$gene <- NA
    removeList <- vector()
    for ( i in seq( 1 , nrow( regions ) ) )
    {
      if ( regions[i,"chromosome_name"] %in% names( xAxisInfo$cLengths ) )
      {
        regions[i,"xmin"] <- regions[i,"start_position"] + xAxisInfo$cLengths[regions[i,"chromosome_name"]]
        regions[i,"xmax"] <- regions[i,"end_position"] + xAxisInfo$cLengths[regions[i,"chromosome_name"]]
        regions[i,"xmean"] <- mean( c( regions[i,"xmin"] , regions[i,"xmax"] ) )
        regions[i,"gene"] <- regions[i,"external_gene_id"]
      }
      else
      {
        removeList <- c( removeList , i )
      }
    }
    if ( length( removeList ) > 0 )
    {
      regions <- regions[-removeList,]
    }
    hexPVals <- hexPVals + geom_rect( data = regions , aes( NULL , NULL , xmin = xmin , xmax = xmax , ymin = 0 , ymax = ymax ) , colour = "red" , alpha = 0.3 )
    hexPVals <- hexPVals + geom_text( data = regions , aes( label = gene , x = xmean , y = ymax ) , hjust = 0 , vjust = -1 , size = 2 , angle = 45 )
  }
  if ( rawChromosomes )
  {
    bleacher <- function( x , bleach = 0 )
    {
      return( ( x * ( 1 - bleach ) ) + bleach )
    }
    chrom.width <- 0.5
    require( grid )
    require( lodplot )
    data( chrom.bands , "chrom.bands" , package = "lodplot" )
    units <- "bases"
    ymin <- 0.975 * max( x$pValue )
    ymax <- 0.995 * max( x$pValue )
    for ( chrom in names( xAxisInfo$xLabelLocs ) )
    {
      chrom <- switch( as.character( chrom ) , `98` = "X" , `99` = "Y" , as.character( chrom ) )
      chromdata <- subset( chrom.bands , chrom.bands$chr == chrom )
      if ( nrow( chromdata ) > 0 )
      {
        lc <- nchar( chromdata$band )
        sel <- !( substr( chromdata$band , lc , lc ) %in% letters )
        chromdata <- chromdata[sel,]
        rm( lc, sel )
        bandpos <- switch( units ,
            cM = chromdata[,c( "cM.top" , "cM.bot" )] ,
            bases = chromdata[,c( "bases.top" , "bases.bot" )] , 
            ISCN = chromdata[,c( "ISCN.top ", "ISCN.bot" )]
        )
        type.b <- match( chromdata$stain , c( "acen" , "gneg" , "gpos" , "gvar" , "stalk" ) )
        bandcol <- gray( bleacher( c( 0.5 , 1 , 0.2 , 0.6 , 0.75 ) ) )[type.b]
        banddens <- c( 30 , -1 , -1 , -1 , 10 )[type.b]
        bandbord <- gray( bleacher( c( 0 , 0 , 0 , 0 , 1 ) ) )[type.b]
        n <- nrow( chromdata )
        centromere <- which( chromdata$arm[-n] != chromdata$arm[-1] )
        idx <- c( 2:( centromere - 1 ) , ( centromere + 2 ):( n - 1) )
        pos.bottom <- 0
        pos.top <- chrom.width
        pos.chrom <- ( 1 + chrom.width ) / 2
        pos.band <- chrom.width + 0.1 * ( 1 - chrom.width )
        bleachblack <- gray( bleacher( 0 ) )
        i <- which( chrom %in% names( xAxisInfo$xLabelLocs ) )
        bandpos <- bandpos + xAxisInfo$cLengths[i]
        xmin <- bandpos[idx,1]
        xmax <- bandpos[idx,2]
        xtmp <- data.frame( xmin = xmin , xmax = xmax , ymin = ymin , ymax = ymax , fill = bandcol[idx] , colour = bandbord[idx] )
        j <- 1
#        for ( j in seq( 1 , length( xmin ) ) )
#        {
        hexChrom <- ggplot( xtmp ) + geom_rect( aes( xmin = xmin , xmax = xmax , ymin = ymin , ymax = ymax , fill = fill , colour = colour ) ) + scale_color_discrete( guide = FALSE ) + scale_fill_discrete( guide = FALSE ) + opts( panel.grid.major = theme_blank() , panel.grid.minor = theme_blank() , panel.background = theme_blank() , panel.axis.x = theme_blank() , panel.axis.y = theme_blank() , axis.ticks.x = theme_blank() , axis.ticks.y = theme_blank() , axis.title.x = theme_blank() , axis.text.x = theme_blank() , axis.title.y = theme_blank() , axis.text.y = theme_blank() )
#        }
#        fg <- frameGrob()
#        fg <- packGrob( fg , grid.rect( x = unit( bandpos[idx,1] , "native" ) , y = unit( pos.bottom , "native" ) , width = unit( bandpos[idx,2] - bandpos[idx, 1] , "native" ) , height = unit( chrom.width , "native" ) , just = c( "left", "bottom" ) , default.units = "native" , gp = gpar( fill = bandcol[idx] , col = bandbord[idx] ) ) )
#        fg <- packGrob( fg , qsGridSemicircle( bandpos[1,2] , pos.bottom , chrom.width , bandpos[1,2] - bandpos[1,1] , default.units = "native" , gp = gpar( fill = bandcol[1] , col = bandbord[1] ) ) )
#        fg <- packGrob( fg , quantsmooth:::qs.grid.semicircle( bandpos[n,1] , pos.bottom , chrom.width , bandpos[n,2] - bandpos[n,1] , default.units = "native" , gp = gpar( fill = bandcol[n] , col = bandbord[n] ) ) )
#        fg <- packGrob( fg , quantsmooth:::qs.grid.semicircle( bandpos[centromere,1] , pos.bottom , chrom.width , bandpos[centromere,2] - bandpos[centromere,1] , default.units = "native" , gp = gpar( fill = bandcol[centromere], col = bandbord[centromere] ) ) )
#        fg <- packGrob( fg , quantsmooth:::qs.grid.semicircle( bandpos[centromere + 1,2] , pos.bottom , chrom.width , bandpos[centromere + 1,2] - bandpos[centromere + 1,1] , default.units = "native" , gp = gpar( fill = bandcol[centromere +1] , col = bandbord[centromere + 1] ) ) )
#        fg <- packGrob( fg , grid.circle( bandpos[centromere,2] , pos.bottom + chrom.width / 2 , unit( chrom.width * 0.3 , "npc" ) , default.units = "native" , gp = gpar( col = bleachblack , fill = "white" , lwd = 2 ) ) )
#        hexPVals <- hexPVals + annotation_custom( grob = fg , xmin = xAxisInfo$cLengths[i] , xmax = xAxisInfo$cLengths[i + 1] , ymin = ymin , ymax = ymax )
        hexPVals <- hexPVals + annotation_custom( grob = ggplotGrob( hexChrom ) , xmin = xAxisInfo$cLengths[i] , xmax = xAxisInfo$cLengths[i + 1] , ymin = ymin , ymax = ymax )
      }
    }
  }
  print( hexPVals )
  invisible( hexPVals )
}
