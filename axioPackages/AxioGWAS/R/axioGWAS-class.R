#glm:
# [1] "coefficients"      "residuals"         "fitted.values"     "effects"           "R"                 "rank"
# [7] "qr"                "family"            "linear.predictors" "deviance"          "aic"               "null.deviance"
#[13] "iter"              "weights"           "prior.weights"     "df.residual"       "df.null"           "y"
#[19] "converged"         "boundary"          "model"             "call"              "formula"           "terms"
#[25] "data"              "offset"            "control"           "method"            "contrasts"         "xlevels"
#
#lm:
#[1] "coefficients"  "residuals"     "effects"       "rank"          "fitted.values" "assign"        "qr"            "df.residual"
#[9] "contrasts"     "xlevels"       "call"          "terms"         "model"
#
#lme:
# [1] "modelStruct"  "dims"         "contrasts"    "coefficients" "varFix"       "sigma"        "apVar"        "logLik"
# [9] "numIter"      "groups"       "call"         "terms"        "method"       "fitted"       "residuals"    "fixDF"
#[17] "na.action"    "data"
#
#multinom:
# [1] "n"             "nunits"        "nconn"         "conn"          "nsunits"       "decay"         "entropy"       "softmax"
# [9] "censored"      "value"         "wts"           "convergence"   "fitted.values" "residuals"     "lev"           "call"
#[17] "terms"         "weights"       "deviance"      "rank"          "coefnames"     "vcoefnames"    "contrasts"     "xlevels"
#[25] "edf"           "AIC"


#' axioGWAS class
#'
#' The Class of an object to store all the data and parameters necessary to run a GWAS, as well as present the results. To create the object, use the new('axioGWAS') command.
#'
#' @examples
#' \dontrun{
#' examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
#' exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B","A/A","A/B","B/B","A/B","A/A","B/B",
#' "A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B","A/A","A/B","B/B","A/B","A/A","B/B",
#' "A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B","A/A","A/B","A/A","A/B","B/B","A/B","A/A","B/B"),
#' nrow=9,byrow=TRUE,dimnames=list(c("snp1","snp2","snp3","snp4","snp5","snp6","snp7","snp8","snp9")))
#' exampleSNPInfo<-data.frame(Chr=c(1,1,1,2,2,2,3,3,3),Pos=c(1,2,3,1,2,3,1,2,3),row.names = c("snp1","snp2","snp3","snp4","snp5","snp6","snp7","snp8","snp9"))
#' exampleModel<-formula("var1~var2+var3+gender")
#'
#' axioGWAS<-new('axioGWAS')
#' axioGWAS$genotypes<-exampleSNP
#' axioGWAS$phenotypes<-examplePGx
#' axioGWAS$genderVar<-'gender'
#' axioGWAS$interactionTerm<-'var2'
#' axioGWAS$snpFeatures<-exampleSNPInfo
#' axioGWAS$chromosomeVar<-'Chr'
#' axioGWAS$positionVar<-'Pos'
#' axioGWAS$model<-exampleModel
#' axioGWAS$database_path<-"runGWAS.DB"
#' axioGWAS$database_table<-"runGWAS_Test"
#' axioGWAS$overwritepreviousEntry<-TRUE
#' axioGWAS$functionArgs<-list(snpInfo=exampleSNPInfo)
#'
#' axioGWAS$runFilter()
#'
#' axioGWAS$runGWAS()
#'
#' gwasRes<-runGWAS(pgx=examplePGx,snp=exampleSNP,Model=exampleModel,func = gwasFunc,interactionTerm="var2",genderVar="gender",databasePath="runGWAS.DB", dbTable="runGWAS_Test2", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
#' gwasRes$adjPVals$SnPName<-rownames(gwasRes$adjPVals)
#' colnames(gwasRes$adjPVals)<-c('Pvalue','AdjustedPvalue','SNPName')
#' gwasRes$adjPVals<-gwasRes$adjPVals[,c('SNPName','Pvalue','AdjustedPvalue')]
#' identical(axioGWAS$pvalues,gwasRes$adjPVals)
#'
#' axioGWAS$manhattanPlot(hexBins = FALSE)
#'
#' snp1Model<-axioGWAS$extractSNPResult("snp1")
#'
#'
#' }
#' @import AxioSerializer
#' @import ggplot2
#' @import stats
#' @import snowfall
#' @import ff
#' @import dbdf
#' 
#' @export

setRefClass(
  #class name
  Class = "axioGWAS",
  
  #add fields/slots
  fields = list(
    #basic required inputs
    genotypes = "ANY",
    phenotypes = "ANY",
    snpFeatures = "ANY",
    chromosomeVar = "character",
    positionVar = "character",
    significantCutOff = "numeric",
    
    #basic outputs
    significantSNPs = "character",
    pvalues = "data.frame",
    
    #inputs for running GWAS
    func = "function",
    functionArgs = "list",
    model = "formula",
    database_path = "character",
    database_table = "character",
    overwritepreviousEntry = "logical",
    genderVar = "character",
    interactionTerm = "character",
    adjust = "character",
    
    #Snowfall settings
    snowFall = "logical",
    ncores = "numeric",
    tempFile = "character",
    
    #memory Optimization Scheme
    memOptScheme = "character"
  ),
  
  methods = list(
    #Set initial values for some of the fields
    initialize = function(...,
                          chromosomeVar = "Chromosome" ,
                          positionVar = "Position" ,
                          significantCutOff = 0.05,
                          func = AxioGWAS::gwasFunc,
                          overwritepreviousEntry = FALSE,
                          adjust = "BY",
                          snowFall = FALSE,
                          memOptScheme = "none",
                          tempFile = tempfile()) {
      chromosomeVar <<- chromosomeVar
      positionVar <<- positionVar
      significantCutOff <<- significantCutOff
      func <<- func
      overwritepreviousEntry <<- overwritepreviousEntry
      adjust <<- adjust
      snowFall <<- snowFall
      memOptScheme <<- match.arg(memOptScheme,c("ff","dbDF","none"))
      tempFile <<- tempFile
    },
    
    #######################################################################################################
    #### Functions to convert what is beleived to be the slowest part/memory eating part of the gwas   #### 
    #### into a pointer                                                                                #### 
    #######################################################################################################
    optimizeMemory = function(...){
      
      args<-list(...)

      method<-.self$memOptScheme
      if("method"%in%names(args)){
        method=args$method
      }
      
      verbose<-FALSE
      if("verbose"%in%names(args)){
        verbose=args$verbose
      }
      
      method<-match.arg(method,c("ff","dbDF","none"))
      #Convert genotypes
      if(inherits(genotypes,c("data.frame","matrix","data.table"))){
        if(method=="ff"){
          genotypes<<-as.ffdf(as.data.frame(genotypes))
        }else if(method=="dbDF"){
          genotypes<<-as.dbDF(as.data.frame(genotypes))
        }
      }else if(is.ffdf(genotypes)){
        if(method=="ff"){
          if(verbose) message("Genotype object is already an ffdf")
        }else if(method=="dbDF"){
          genotypes<<-as.dbDF(as.data.frame(genotypes))
        }
      }else if(is.dbDF(genotypes)){
        if(method=="ff"){
          genotypes<<-as.ffdf(genotypes[,])
        }else if(method=="dbDF"){
          if(verbose) message("Genotype object is already an dbDF")
        }
      }else{
        warning(paste("Did not recognize genotype object type, could not not optimize storage:",class(genotypes)[1]))
      }
      
      #Convert phenotypes
      if(inherits(phenotypes,c("data.frame","matrix","data.table"))){
        if(method=="ff"){
          phenotypes<<-as.ffdf(as.data.frame(phenotypes))
        }else if(method=="dbDF"){
          phenotypes<<-as.dbDF(as.data.frame(phenotypes))
        }
      }else if(is.ffdf(phenotypes)){
        if(method=="ff"){
          if(verbose) message("phenotypes object is already an ffdf")
        }else if(method=="dbDF"){
          phenotypes<<-as.dbDF(as.data.frame(phenotypes))
        }
      }else if(is.dbDF(phenotypes)){
        if(method=="ff"){
          phenotypes<<-as.ffdf(phenotypes[,])
        }else if(method=="dbDF"){
          if(verbose) message("phenotypes object is already an dbDF")
        }
      }else{
        warning(paste("Did not recognize phenotypes object type, could not not optimize storage:",class(phenotypes)[1]))
      }
      
      #Convert snpFeat
      if(inherits(snpFeatures,c("data.frame","matrix","data.table"))){
        if(method=="ff"){
          snpFeatures<<-as.ffdf(as.data.frame(snpFeatures))
        }else if(method=="dbDF"){
          snpFeatures<<-as.dbDF(as.data.frame(snpFeatures))
        }
      }else if(is.ffdf(snpFeatures)){
        if(method=="ff"){
          if(verbose) message("snpFeatures object is already an ffdf")
        }else if(method=="dbDF"){
          snpFeatures<<-as.dbDF(as.data.frame(snpFeatures))
        }
      }else if(is.dbDF(snpFeatures)){
        if(method=="ff"){
          snpFeatures<<-as.ffdf(snpFeatures[,])
        }else if(method=="dbDF"){
          if(verbose) message("snpFeatures object is already an dbDF")
        }
      }else{
        warning(paste("Did not recognize snpFeatures object type, could not not optimize storage:",class(snpFeatures)[1]))
      }
      
    },
    
    ##########################################################################################
    #### Filter SNPS on call rate, minor allele frequency, and hardy weinberg equilibrium ####
    ##########################################################################################
    runFilter = function(crCutOff = 0.95 ,
                         mafCutOff = 0.05 ,
                         hweCutOff = 0.05 ,
                         chromosomeVar = NULL ,
                         adjust = NULL,
                         returnFilter = FALSE ,
                         returnAll = FALSE) {
      if (snowFall) {
        sfInit(parallel = TRUE, cpus = ncores)
        on.exit(sfStop()) # if there is an erorr, it will automatically end the cluster
        sfLibrary(AxioGWAS)
      }
      
      SNPFilter <- AxioGWAS::runsnpFilter(
        x = genotypes,
        crCutOff = crCutOff,
        mafCutOff = mafCutOff,
        hweCutOff = hweCutOff,
        Chr = setNames(as.character(snpFeatures[rownames(genotypes), ifelse(!is.null(chromosomeVar) ,
                                                                            chromosomeVar ,
                                                                            .self$chromosomeVar)]), rownames(genotypes)),
        gender = as.character(phenotypes[colnames(genotypes), genderVar]),
        adjust = ifelse(!is.null(adjust) , adjust , .self$adjust),
        returnAll = returnAll,
        snowFall = snowFall
      )
      
      if (!returnFilter && !returnAll) {
        originalDim<-dim(genotypes)
        genotypes <<- genotypes[SNPFilter,]
        newDim<-dim(genotypes)
        message(paste("Filtering removed",originalDim[1]-newDim[1],"SNPs, the number of SNPs ready to be analyzed is:", newDim[1]))
        
      } else{
        
        return(SNPFilter)
        
      }
      
      
    },
    
    
    #############################################################################
    #### Functions to Create formula to be used, return the base formula, etc ####
    #############################################################################
    returnModel = function() {
      if (length(.self$model) == 0) {
        stop("No model provided")
      } else{
        message(paste(
          "Model used for the GWAS is as follows:",
          deparse(.self$model),
          "",
          sep = " "
        ))
      }
    },
    
    generateModel = function(depVar,
                             indepVar = NULL,
                             pCutOff = 0.1) {
      if(dim(phenotypes)[1]==0){
        error('The phenotypes field is empty. Please have a valid data.frame entry for the phenotypes field before attempting to generate a model.')
      }
      
      if (is.null(indepVar)) {
        indepVar <- setdiff(colnames(phenotypes), depVar)
      }
      tmpModel <- try(findModel(phenotypes, depVar, indepVar, pCutOff))
      if (!is(tmpModel, "try-error")) {
        model <<- formula(tmpModel)
        returnModel()
        
      } else{
        warning("Unable to generate a model from supplied depVar and indepVar values")
      }
      
    },
    
    
    ##########################################################################
    #### Functions to run GWAS, and to filter results for significant SNPS ####
    ##########################################################################
    
    runGWAS = function(optimizeMemory=TRUE) {
      
      if(optimizeMemory){
        .self$optimizeMemory()
      }
      
      
      if(file.exists(paste0(tempFile,".log"))){
        file.remove(paste0(tempFile,".log"))
      }
     
      if (snowFall) {
        sfInit(parallel = TRUE, cpus = ncores)
        on.exit(sfStop())
        sfLibrary(AxioGWAS)
        
        dbConnection <- dbConnect(SQLite(), database_path)
        on.exit({if(dbIsValid(dbConnection))dbDisconnect(dbConnection)},add=TRUE)
        
        if (!(database_table %in% dbListTables(dbConnection))) {
          dbExecute(
              dbConnection,
              paste(
                  "CREATE TABLE",
                  database_table,
                  "(OBJECTNAMES varchar PRIMARY KEY,RAWDATA varchar)"
              )
          )
          dbExecute(
              dbConnection,
              paste(
                  "CREATE TABLE",
                  paste( "Anova" , database_table , sep = "_" ),
                  "(OBJECTNAMES varchar PRIMARY KEY,RAWDATA varchar)"
              )
          )
        }
        dbDisconnect(dbConnection)
        
      }
      
      runGWASCommand <-
        'gwasRES <- AxioGWAS::runGWAS( pgx = phenotypes , snp = genotypes , Model = model , func = func , interactionTerm = interactionTerm , adjust = adjust , snowFall = snowFall , genderVar = genderVar , dbPath = database_path , dbTABLENAME = database_table , overwriteEntry = overwritepreviousEntry , tempFile = tempFile'
      if (length(functionArgs) > 0)
      {
        unnecessaryArgs <- setdiff( names( functionArgs ) , names( formals( func ) ) )
        
        if ( length( unnecessaryArgs ) > 0 )
        {
          if ( any( unnecessaryArgs %in% "correlationForm" ) )
          {
            if ( any( names( formals( func ) ) %in% c( "corrFun" , "phiVal" , "corrForm" ) ) )
            {
              unnecessaryArgs <- unnecessaryArgs[-which( unnecessaryArgs == "correlationForm" )]
            }
          }
          if ( length( unnecessaryArgs ) > 0 )
          {
            functionArgs <- functionArgs[-which( names( functionArgs ) %in% unnecessaryArgs )]
          }
        }
        
        additionalArgs <- c()
        for ( i in seq_along( names( functionArgs ) ) )
        {
          if ( inherits( functionArgs[[i]] , "corStruct" ) )
          {
            funArgList <- getCorList( functionArgs[[i]] )
            for ( j in seq_along( names( funArgList ) ) )
            {
              additionalArgs <- c( additionalArgs , paste( names( funArgList )[j] , ifelse( inherits( funArgList[[j]] , c( "formula" , "call" ) ) , paste( "\"" , deparse( funArgList[[j]] ) , "\"" , sep = "" ) , paste( "\"" , funArgList[[j]] , "\"" , sep = "" ) ) , sep = " = "  ) )
            }
          }
          else if( inherits( functionArgs[[i]] , "glmerControl" ) ){
            
          }
          else if(inherits( functionArgs[[i]] , c("data.frame","ff","dbDF","matrix","data.table" ) ) ){
            additionalArgs <- c( additionalArgs , paste( names( functionArgs) [i] , paste0("functionArgs[[",i,"]]") , sep = " = "  ) )
          }
          else
          {
            additionalArgs <- c( additionalArgs , paste( names( functionArgs )[i] , ifelse( inherits( functionArgs[[i]] , c( "formula" , "call" ) ) , paste( "\"" , deparse( functionArgs[[i]] ) , "\"" , sep = "" ) , paste( "\"" , functionArgs[[i]] , "\"" , sep = "" ) ) , sep = " = "  ) )
          }
        }
        runGWASCommand <- paste0( c( runGWASCommand , additionalArgs ) , collapse = " , " )
      }
      runGWASCommand <- paste0( runGWASCommand , ")" )
      
      message(paste(
        "Starting GWAS with command:",
        gsub("gwasRES <-", "", runGWASCommand)
      ))
      
      eval(parse(text = runGWASCommand))
      
      colnames(gwasRES)<-c('Pvalue','SNPName','AdjustedPvalue')

      pvalues <<- gwasRES[,c('SNPName','Pvalue','AdjustedPvalue')]

      findSignificantSNPS()
      
      if (snowFall) {
        sfStop()
      }
      invisible( gwasRES )
    },
    
    #find significant snps from GWAS
    findSignificantSNPS = function(cutoff = NULL) {
      if(is.null(cutoff)){
        cutoff<-significantCutOff
      }
      if(length(which(pvalues$AdjustedPvalue < cutoff))>0){
      significantSNPs <<-
        rownames(pvalues)[which(pvalues$AdjustedPvalue < cutoff)]
      }
    },
    
    ####################################################
    #### Functions to generate Plots of GWAS results ####
    ####################################################
    
    manhattanPlot = function(snpFeat = NULL ,
        chromVar = NULL,
        posVar = NULL,
        chromosome = NULL,
        plotSignificantSNPS = TRUE,
        hexBins = TRUE ) {
      if (is.null(snpFeat)) {
        snpFeat <- .self$snpFeatures
        if (is.null(chromVar)) {
          chromVar <- .self$chromosomeVar
        }
        if (is.null(posVar)) {
          posVar <- .self$positionVar
        }
      }
      
      ManhattanInfo <-
          data.frame(
              P.Value = as.numeric(as.character(pvalues$Pvalue)),
              Position = as.numeric(as.character(snpFeat[rownames(pvalues), posVar])),
              Chromosome = as.character(snpFeat[rownames(pvalues), chromVar]),
              AdjPValue = as.numeric(as.character(pvalues$AdjustedPvalue)),
              LogP = -log10(as.numeric(as.character(pvalues$Pvalue))),
              LogAdjpvalue = -log10(as.numeric(as.character(pvalues$AdjustedPvalue)))
          )
      rownames(ManhattanInfo) <- rownames(pvalues)
      
      if (is.null(chromosome)) {
        chromosome <- unique(as.character(ManhattanInfo[, "Chromosome"]))
        suppressWarnings(numericChromosomes <-
                chromosome[which(!is.na(as.numeric(chromosome)))])
        numericChromosomes <-
            numericChromosomes[order(as.numeric(numericChromosomes))]
        suppressWarnings(characterChromosomes <-
                chromosome[which(is.na(as.numeric(chromosome)))])
        chromosome <- c(numericChromosomes, characterChromosomes)
      }
      
      ChromosomePositions <- list()
      ChromosomePositions[[chromosome[1]]] <- c(Start = 0,
          LabelPosition = max(ManhattanInfo[which(ManhattanInfo$Chromosome == chromosome[1]), 'Position'],na.rm = TRUE)/2)
      if (length(chromosome) > 1) {
        for (i in seq(2,length(chromosome))) {
          Chrom <- chromosome[i]
          prevChrom <- chromosome[i - 1]
          ChromosomePositions[[Chrom]] <-
              c(Start = max(ManhattanInfo[which(ManhattanInfo$Chromosome == prevChrom), 'Position']),
                  LabelPosition = max(ManhattanInfo[which(ManhattanInfo$Chromosome == Chrom), 'Position'],na.rm = TRUE)/2)
        }
      }
      ChromosomePositions <- do.call('rbind', ChromosomePositions)
      ChromosomePositions[, 'Start'] <-
          cumsum(ChromosomePositions[, 'Start'])
      
      ChromosomeLabs<-setNames(ChromosomePositions[, 'Start']+ChromosomePositions[,'LabelPosition'],
          rownames(ChromosomePositions))
      ChromosomePositions <-
          setNames(ChromosomePositions[, 'Start'], rownames(ChromosomePositions))
      
      ManhattanInfo <-
          ManhattanInfo[which(ManhattanInfo$Chromosome %in% names(ChromosomePositions)), ]
      ManhattanInfo$position <-
          ManhattanInfo$Position + ChromosomePositions[as.character(ManhattanInfo$Chromosome)]
      if (any(is.na(ManhattanInfo$P.Value))) {
        ManhattanInfo <- ManhattanInfo[!is.na(ManhattanInfo$"P.Value"), ]
      }
      ManhattanInfo$Group <- NA
      ManhattanInfo$Group[which(as.numeric(ManhattanInfo$Chromosome) %% 2 == 0)] <-
          "red"
      ManhattanInfo$Group[which(as.numeric(ManhattanInfo$Chromosome) %% 2 != 0)] <-
          "black"
      
      ManhattanInfo$Chromosome <-
          factor(ManhattanInfo$Chromosome, levels = chromosome)
      
      
      
      ggManhattan <-
          ggplot(ManhattanInfo) + geom_point(aes(
                  x = position,
                  y = LogAdjpvalue,
                  size = 3.5,
                  alpha = 1 / 3,
                  colour = Group
              ),
              show.legend = FALSE)
      ggManhattan <- ggManhattan + theme_bw(base_size = 15)
      ggManhattan <-
          ggManhattan + scale_x_continuous(breaks = as.numeric(ChromosomeLabs),
              labels = names(ChromosomeLabs) ) + geom_vline(xintercept = as.numeric(ChromosomePositions),
              color = "grey")
      ggManhattan <-
          ggManhattan + geom_hline(
              yintercept = -log10(.05),
              linetype = 1,
              col = 'red',
              lwd = 1.5
          ) + geom_text(aes(
                  label = "p-Value = 0.05" ,
                  x = 0 ,
                  y = 2.2 ,
                  hjust = 0,
                  size = 5
              ),
              show.legend = FALSE)
      ggManhattan <-
          ggManhattan + geom_segment(aes(
                  x = 0 ,
                  y = 2 ,
                  xend = 0 ,
                  yend = (-log10(.05) + .01)
              ) , arrow = arrow(length = unit(0.25 , "cm")))
      ggManhattan <-
          ggManhattan + ggtitle('Manhattan Plot') + xlab("Chromosomal Position") + ylab('-log10(adj p-Value)')
      
      
      
      if (plotSignificantSNPS) {
        if (length(which(significantSNPs%in%rownames(ManhattanInfo))) > 0) {
          IntSnPManhat <- ManhattanInfo[significantSNPs, ]
          if (any(is.na(IntSnPManhat$P.Value))) {
            IntSnPManhat <- IntSnPManhat[!is.na(IntSnPManhat$"P.Value"), ]
          }
          ggManhattan <-
              ggManhattan + geom_point(
                  data = IntSnPManhat,
                  aes(
                      x = position,
                      y = LogAdjpvalue,
                      color = "red",
                      size = 10
                  ),
                  alpha = 1 / 3,
                  show.legend = FALSE
              )
        }
        
      }
      
      if(hexBins){
        #removes the points below the significant value, and recreates them as hexbins
        if(length(which(ggManhattan$data$AdjPValue<=0.05))>0){
          ggManhattan$data<-ggManhattan$data[which(ggManhattan$data$AdjPValue<0.05),]
        }else{
          ggManhattan$layers[[1]]<-NULL
        }
        ggManhattan<-ggManhattan+geom_hex(data=ManhattanInfo[which(ManhattanInfo$AdjPValue>0.05),],aes(x = position,y = LogAdjpvalue))
        
      }
      
      ggManhattan + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )
    },
    
    ########################################################################
    #### Functions to extract individual linear models from GWAS results ####
    ########################################################################
    
    extractSNPResult = function( snpName , verbose = TRUE )
    {
      if ( verbose )
      {
        message( paste( "Extracting model for" , snpName ) )
      }
      dbConnection <- dbConnect( SQLite() , database_path )
      listModels <- dbGetQuery( dbConnection , paste( "SELECT OBJECTNAMES FROM" , database_table ) )$OBJECTNAMES
      if ( length( grep( snpName , listModels ) ) == 0 )
      {
        snpName <- as.character( which( rownames( genotypes ) == snpName ) )
      }
      if ( length( grep( snpName , listModels ) ) != 0 )
      {
        modeltoExtract <- listModels[grep( snpName , listModels )]
        dbDisconnect( dbConnection )
        return( readObjectFromTable( modeltoExtract , database_table , database_path ) )
      }
      else
      {
        stop( paste( "The regression data for" , snpName , "could not be found. Check spelling or if SNP passed filtering" ) )
      }
    },
    
    ########################################################################
    #### Functions to extract individual ANOVA tables from GWAS results ####
    ########################################################################
    
    extractAnova = function( snpName , verbose = TRUE )
    {
      if ( verbose )
      {
        message( paste( "Extracting ANOVA table for" , snpName ) )
      }
      dbConnection <- dbConnect( SQLite() , database_path )
      listModels <- dbGetQuery( dbConnection , paste( "SELECT OBJECTNAMES FROM" , paste( "Anova" , database_table , sep = "_" ) ) )$OBJECTNAMES
      if ( length( grep( snpName , listModels ) ) == 0 )
      {
        snpName <- as.character( which( rownames( genotypes ) == snpName ) )
      }
      if ( length( grep( snpName , listModels ) ) != 0 )
      {
        modeltoExtract <- listModels[grep( snpName , listModels )]
        dbDisconnect( dbConnection )
        return( readObjectFromTable( modeltoExtract , paste( "Anova" , database_table , sep = "_" ) , database_path ) )
      }
      else
      {
        stop( paste( "The ANOVA table for" , snpName , "could not be found. Check spelling or if SNP passed filtering" ) )
      }
    },
    
    ########################################################################
    #### Functions to run Anova from GWAS results ####
    ########################################################################
    
    runAnova = function( snpName , ... )
    {
      message( paste( "Running Anova for" , snpName ) )
      snpOut <- extractSNPResult( snpName , FALSE )
      if ( inherits( snpOut , c( "coxph" , "lme" ) ) )
      {
        return( anova( snpOut , ... ) )
      }
      else
      {
        return( Anova( snpOut , ... ) )
      }
    }
  )
)
