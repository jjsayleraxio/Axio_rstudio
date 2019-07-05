runControlReport <- function( externalData = NULL , SQL = NULL , reportFile = NULL , reportPath = NULL , fileList = NULL , failList = NULL , dna10ThreshLow = NULL , dna100ThreshLow = NULL , rna10ThreshLow = NULL , rna100ThreshLow = NULL , dna10ThreshHigh = NULL , dna100ThreshHigh = NULL , rna10ThreshHigh = NULL , rna100ThreshHigh = NULL , TRIM = TRUE , dnaLCL = 0.025 , dnaUCL = 0.975 , rnaLCL = 0.025 , rnaUCL = 0.975 , dateMin = NULL , dateMax = NULL , excludeList = NULL , annotationList = NULL )
{
  if ( !is.null( fileList ) )
  {
    controlData <- readControls( fileList )
  }

  if ( !is.null( externalData ) )
  {
    if ( !is.null( fileList ) )
    {
      for ( i in names( externalData ) )
      {
        controlData[[i]] <- rbind( controlData[[i]] , externalData[[i]] )
      }
    }
    else
    {
      controlData <- externalData
    }
  }

  if ( !is.null( SQL ) )
  {
    if ( !is.null( fileList ) )
    {
      for ( i in names( externalData ) )
      {
        controlData[[i]] <- rbind( controlData[[i]] , externalData[[i]] )
      }
    }
    else
    {
      controlData <- externalData
    }
  }

  if ( !is.null( failList ) )
  {
    failData <- readFails( Files = failList )
    controlData[["DNA10"]] <- rbind( controlData[["DNA10"]] , failData[["DNA10"]] )
    controlData[["DNA100"]] <- rbind( controlData[["DNA100"]] , failData[["DNA100"]] )
    controlData[["RNA10"]] <- rbind( controlData[["RNA10"]] , failData[["RNA10"]] )
    controlData[["RNA100"]] <- rbind( controlData[["RNA100"]] , failData[["RNA100"]] )
  }

  if ( !is.null( dateMin ) )
  {
    controlData[["DNA10"]] <- controlData[["DNA10"]][which( controlData[["DNA10"]]$"QC_Date" > dateMin ),]
    controlData[["DNA100"]] <- controlData[["DNA100"]][which( controlData[["DNA100"]]$"QC_Date" > dateMin ),]
    controlData[["RNA10"]] <- controlData[["RNA10"]][which( controlData[["RNA10"]]$"QC_Date" > dateMin ),]
    controlData[["RNA100"]] <- controlData[["RNA100"]][which( controlData[["RNA100"]]$"QC_Date" > dateMin ),]
  }
  if ( !is.null( dateMax ) )
  {
    controlData[["DNA10"]] <- controlData[["DNA10"]][which( controlData[["DNA10"]]$"QC_Date" < dateMax ),]
    controlData[["DNA100"]] <- controlData[["DNA100"]][which( controlData[["DNA100"]]$"QC_Date" < dateMax ),]
    controlData[["RNA10"]] <- controlData[["RNA10"]][which( controlData[["RNA10"]]$"QC_Date" < dateMax ),]
    controlData[["RNA100"]] <- controlData[["RNA100"]][which( controlData[["RNA100"]]$"QC_Date" < dateMax ),]
  }

  if ( !is.null( excludeList ) )
  {
    controlData <- trimDates( controlData , excludeList )
  }
  
  if ( !is.null( dna10ThreshLow ) )
  {
    controlData[["DNA10"]] <- controlData[["DNA10"]][which( controlData[["DNA10"]]$"10ng" > dna10ThreshLow ),]
  }
  
  if ( !is.null( dna10ThreshHigh ) )
  {
    controlData[["DNA10"]] <- controlData[["DNA10"]][which( controlData[["DNA10"]]$"10ng" < dna10ThreshLow ),]
  }
  
  if ( !is.null( dna100ThreshLow ) )
  {
    controlData[["DNA100"]] <- controlData[["DNA100"]][which( controlData[["DNA100"]]$"100ng" > dna100ThreshLow ),]
  }
  
  if ( !is.null( dna100ThreshHigh ) )
  {
    controlData[["DNA100"]] <- controlData[["DNA100"]][which( controlData[["DNA100"]]$"100ng" < dna100ThreshHigh ),]
  }
  
  if ( !is.null( rna10ThreshLow ) )
  {
    controlData[["RNA10"]] <- controlData[["RNA10"]][which( controlData[["RNA10"]]$"10ng" > rna10ThreshLow ),]
  }
  
  if ( !is.null( rna10ThreshHigh ) )
  {
    controlData[["RNA10"]] <- controlData[["RNA10"]][which( controlData[["RNA10"]]$"10ng" < rna10ThreshHigh ),]
  }
  
  if ( !is.null( rna100ThreshLow ) )
  {
    controlData[["RNA100"]] <- controlData[["RNA100"]][which( controlData[["RNA100"]]$"100ng" > rna100ThreshLow ),]
  }
  
  if ( !is.null( rna100ThreshHigh ) )
  {
    controlData[["RNA100"]] <- controlData[["RNA100"]][which( controlData[["RNA100"]]$"100ng" < rna100ThreshHigh ),]
  }

  if ( TRIM )
  {
    dna10lcl <- quantile( controlData[["DNA10"]]$"10ng" , dnaLCL )
    dna10ucl <- quantile( controlData[["DNA10"]]$"10ng" , dnaUCL )
    dna100lcl <- quantile( controlData[["DNA100"]]$"100ng" , dnaLCL )
    dna100ucl <- quantile( controlData[["DNA100"]]$"100ng" , dnaUCL )
    
    rna10lcl <- quantile( controlData[["RNA10"]]$"10ng" , rnaLCL )
    rna10ucl <- quantile( controlData[["RNA10"]]$"10ng" , rnaUCL )
    rna100lcl <- quantile( controlData[["RNA100"]]$"100ng" , rnaLCL )
    rna100ucl <- quantile( controlData[["RNA100"]]$"100ng" , rnaUCL )
    
    controlData[["DNA10"]] <- controlData[["DNA10"]][which( controlData[["DNA10"]]$"10ng" > dna10lcl & controlData[["DNA10"]]$"10ng" < dna10ucl ),]
    controlData[["DNA100"]] <- controlData[["DNA100"]][which( controlData[["DNA100"]]$"100ng" > dna100lcl & controlData[["DNA100"]]$"100ng" < dna100ucl ),]
    controlData[["RNA10"]] <- controlData[["RNA10"]][which( controlData[["RNA10"]]$"10ng" > rna10lcl & controlData[["RNA10"]]$"10ng" < rna10ucl ),]
    controlData[["RNA100"]] <- controlData[["RNA100"]][which( controlData[["RNA100"]]$"100ng" > rna100lcl & controlData[["RNA100"]]$"100ng" < rna100ucl ),]
  }
  
  histRNADataList10 <- runBootDensities( controlData , "RNA" , date = "QC_Date" , level = "10ng" )
  histRNADataList100 <- runBootDensities( controlData , "RNA" , date = "QC_Date" , level = "100ng" )
  histDNADataList10 <- runBootDensities( controlData , "DNA" , date = "QC_Date" , level = "10ng" )
  histDNADataList100 <- runBootDensities( controlData , "DNA" , date = "QC_Date" , level = "100ng" )

  histDenDNA100Plot <- plotDensities( histDNADataList100 , level = "DNA 100ng" , pathRoot = workRoot , main = "Density by Quarter" )
  histDenDNA10Plot <- plotDensities( histDNADataList10 , level = "DNA 10ng" , pathRoot = workRoot , main = "Density by Quarter" )
  histDenRNA100Plot <- plotDensities( histRNADataList100 , level = "RNA 100ng" , pathRoot = workRoot , main = "Density by Quarter" )
  histDenRNA10Plot <- plotDensities( histRNADataList10 , level = "RNA 10ng" , pathRoot = workRoot , main = "Density by Quarter" )

  histDenDNAlim10Plot <- plotLimitWidths( histDNADataList10 , pathRoot = workRoot , level = "10ng" )
  histDenDNAlim100Plot <- plotLimitWidths( histDNADataList100 , pathRoot = workRoot , level = "100ng" )
  histDenRNAlim10Plot <- plotLimitWidths( histRNADataList10 , pathRoot = workRoot , level = "10ng" )
  histDenRNAlim100Plot <- plotLimitWidths( histRNADataList100 , pathRoot = workRoot , level = "100ng" )

  histDenDNACen10Plot <- plotLimitCenters( histDNADataList10 , pathRoot = workRoot , level = "10ng" , PNG = TRUE )
  histDenDNACen100Plot <- plotLimitCenters( histDNADataList100 , pathRoot = workRoot , level = "100ng" , PNG = TRUE )
  histDenRNACen10Plot <- plotLimitCenters( histRNADataList10 , pathRoot = workRoot , level = "10ng" , PNG = TRUE )
  histDenRNACen100Plot <- plotLimitCenters( histRNADataList100 , pathRoot = workRoot , level = "100ng" , PNG = TRUE )

  controlDNABoxPlot <- boxPlotWithSmoothLine( controlData , dnarna = "DNA" , main = "Predicted Concentration versus Date" , lcl10 = dna10lcl , lcl100 = dna100lcl , ucl10 = dna10ucl , ucl100 = dna100ucl , dateMin = dateMin , dateMax = dateMax , annotationList = annotationList[["DNA"]] )
  controlRNABoxPlot <- boxPlotWithSmoothLine( controlData , dnarna = "RNA" , main = "Predicted Concentration versus Date" , lcl10 = rna10lcl , lcl100 = rna100lcl , ucl10 = rna10ucl , ucl100 = rna100ucl , dateMin = dateMin , dateMax = dateMax , annotationList = annotationList[["RNA"]] )

  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesDNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenDNA100Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesDNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenDNA10Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesRNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenRNA100Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesRNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenRNA10Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthDNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenDNAlim10Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthDNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenDNAlim100Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthRNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenRNAlim100Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthRNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenRNAlim10Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterDNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenDNACen100Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterDNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenDNACen10Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterRNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenRNACen100Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterRNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( histDenRNACen10Plot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "controlDNABoxPlot" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( controlDNABoxPlot )
  dev.off()
  png( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "controlRNABoxPlot" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , width = 9.75 , height = 5.8 , units = "in" , res = 125 )
  print( controlRNABoxPlot )
  dev.off()

  controlWorkbook <- loadWorkbook( file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( reportFile , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "xlsx" , sep = "." ) ) , create = TRUE )

  createSheet( controlWorkbook , name = "Limits" )
  writeWorksheet( object = controlWorkbook , data = data.frame( type = "DNA 10 ng/ul" ) , sheet = "Limits" , header = FALSE )
  tRow <- 2
  dna10Limits <- makeLegendFrame( histDNADataList10[["bootPercentiles"]] )
  dna100Limits <- makeLegendFrame( histDNADataList100[["bootPercentiles"]] )
  rna10Limits <- makeLegendFrame( histRNADataList10[["bootPercentiles"]] )
  rna100Limits <- makeLegendFrame( histRNADataList100[["bootPercentiles"]] )
  writeWorksheet( object = controlWorkbook , data = dna10Limits , sheet = "Limits" , startRow = tRow , header = TRUE )
  tRow <- tRow + nrow( dna10Limits ) + 2
  writeWorksheet( object = controlWorkbook , data = data.frame( type = "DNA 100 ng/ul" ) , sheet = "Limits" , startRow = tRow , header = FALSE )
  tRow <- tRow + 1
  writeWorksheet( object = controlWorkbook , data = dna100Limits , sheet = "Limits" , startRow = tRow , header = TRUE )
  tRow <- tRow + nrow( dna100Limits ) + 2
  writeWorksheet( object = controlWorkbook , data = data.frame( type = "RNA 10 ng/ul" ) , sheet = "Limits" , startRow = tRow , header = FALSE )
  tRow <- tRow + 1
  writeWorksheet( object = controlWorkbook , data = rna10Limits , sheet = "Limits" , startRow = tRow , header = TRUE )
  tRow <- tRow + nrow( rna10Limits ) + 2
  writeWorksheet( object = controlWorkbook , data = data.frame( type = "RNA 100 ng/ul" ) , sheet = "Limits" , startRow = tRow , header = FALSE )
  tRow <- tRow + 1
  writeWorksheet( object = controlWorkbook , data = rna100Limits , sheet = "Limits" , startRow = tRow , header = TRUE )

  createSheet( controlWorkbook , name = "DNA_Densities" )
  createName( controlWorkbook , name = "DNA_Densities" , formula = "DNA_Densities!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesDNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "DNA_Densities" , originalSize = TRUE )
  createName( controlWorkbook , name = "DNA_Densities" , formula = "DNA_Densities!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesDNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "DNA_Densities" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "RNA_Densities" )
  createName( controlWorkbook , name = "RNA_Densities" , formula = "RNA_Densities!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesRNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "RNA_Densities" , originalSize = TRUE )
  createName( controlWorkbook , name = "RNA_Densities" , formula = "RNA_Densities!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "densitiesRNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "RNA_Densities" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "DNA_Limit_Widths" )
  createName( controlWorkbook , name = "DNA_Limit_Widths" , formula = "DNA_Limit_Widths!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthDNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "DNA_Limit_Widths" , originalSize = TRUE )
  createName( controlWorkbook , name = "DNA_Limit_Widths" , formula = "DNA_Limit_Widths!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthDNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "DNA_Limit_Widths" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "RNA_Limit_Widths" )
  createName( controlWorkbook , name = "RNA_Limit_Widths" , formula = "RNA_Limit_Widths!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthRNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "RNA_Limit_Widths" , originalSize = TRUE )
  createName( controlWorkbook , name = "RNA_Limit_Widths" , formula = "RNA_Limit_Widths!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitWidthRNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "RNA_Limit_Widths" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "DNA_Limit_Centers" )
  createName( controlWorkbook , name = "DNA_Limit_Centers" , formula = "DNA_Limit_Centers!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterDNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "DNA_Limit_Centers" , originalSize = TRUE )
  createName( controlWorkbook , name = "DNA_Limit_Centers" , formula = "DNA_Limit_Centers!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterDNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "DNA_Limit_Centers" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "RNA_Limit_Centers" )
  createName( controlWorkbook , name = "RNA_Limit_Centers" , formula = "RNA_Limit_Centers!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterRNAControls100" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "RNA_Limit_Centers" , originalSize = TRUE )
  createName( controlWorkbook , name = "RNA_Limit_Centers" , formula = "RNA_Limit_Centers!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "limitCenterRNAControls10" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "RNA_Limit_Centers" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "Box_Plots" )
  writeWorksheet( object = controlWorkbook , data = data.frame( type = "DNA" ) , sheet = "Box_Plots" , startRow = 1 , startCol = 1 , header = FALSE )
  createName( controlWorkbook , name = "Box_Plots" , formula = "Box_Plots!$B$1" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "controlDNABoxPlot" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "Box_Plots" , originalSize = TRUE )
  writeWorksheet( object = controlWorkbook , data = data.frame( type = "RNA" ) , sheet = "Box_Plots" , startRow = 31 , startCol = 1 , header = FALSE )
  createName( controlWorkbook , name = "Box_Plots" , formula = "Box_Plots!$B$31" , overwrite = TRUE )
  addImage( controlWorkbook , filename = file.path( ifelse( is.null( reportPath ) , getwd() , reportPath ) , paste( paste( "controlRNABoxPlot" , format( Sys.Date() , format = "%d%b%Y" ) , sep = "-" ) , "png" , sep = "." ) ) , name = "Box_Plots" , originalSize = TRUE )

  createSheet( controlWorkbook , name = "DNA_100ngul" )
  writeWorksheet( object = controlWorkbook , data = controlData[["DNA100"]] , sheet = "DNA_100ngul" , header = TRUE )
  
  createSheet( controlWorkbook , name = "DNA_10ngul" )
  writeWorksheet( object = controlWorkbook , data = controlData[["DNA10"]] , sheet = "DNA_10ngul" , header = TRUE )
  
  createSheet( controlWorkbook , name = "RNA_100ngul" )
  writeWorksheet( object = controlWorkbook , data = controlData[["RNA100"]] , sheet = "RNA_100ngul" , header = TRUE )
  
  createSheet( controlWorkbook , name = "RNA_10ngul" )
  writeWorksheet( object = controlWorkbook , data = controlData[["RNA10"]] , sheet = "RNA_10ngul" , header = TRUE )

  saveWorkbook( controlWorkbook )
  invisible()
}
