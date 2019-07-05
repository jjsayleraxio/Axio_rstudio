library( XLConnect )

source( "c:/Documents and Settings/dhender1/workspace/BiomekPrep/helperFunctions.R" )

datFile1 <- "C:/Data/BiomekPrep/SC-732418_ribogreentesting_4Oct2011.TXT"
datFile2 <- "C:/Data/BiomekPrep/SC-732417_ribogreentesting_4oct2011.TXT"
datFile3 <- "C:/Data/BiomekPrep/SC-732537_testing_standards_6Oct2011.TXT"
datFile4 <- "C:/Data/BiomekPrep/SC-732538_5oct2011.TXT"
layoutFile <- "C:/Data/BiomekPrep/Ribogreen_Plate_Layout_3oct2011.xls"

datWB <- loadWorkbook( layoutFile )

layout <- readWorksheet( datWB , sheet = "SC-732418_ribogeen_testing" , startRow = 5 , endRow = 21 , startCol = 1 , endCol = 25 , header = TRUE )

Data1 <- read.table( datFile1 , sep = ";" , header = TRUE , row.names = 1 , skip = 2 )
Data1 <- Data1[seq( 1 , 16 , 2 ),]
Data2 <- read.table( datFile2 , sep = ";" , header = TRUE , row.names = 1 , skip = 2 )
Data2 <- Data2[seq( 1 , 16 , 2 ),]
Data3 <- read.table( datFile3 , sep = ";" , header = TRUE , row.names = 1 , skip = 2 )
Data3 <- Data3[seq( 1 , 16 , 2 ),]
Data4 <- read.table( datFile4 , sep = ";" , header = TRUE , row.names = 1 , skip = 2 )
Data4 <- Data4[seq( 1 , 16 , 2 ),]

# Study #732418

layouts732418 <- makeLayouts( layout , Data1 )

ribogreenHigh732418 <- makeRibogreen( layouts732418$oldLayoutHigh , layouts732418$newLayoutHigh , layouts732418$oldDataHigh , layouts732418$newDataHigh )

densityOldHigh732418 <- density( ribogreenHigh732418[["Data"]]$Value[which( ribogreenHigh732418[["Data"]]$set == "Old" )] )
densityNewHigh732418 <- density( ribogreenHigh732418[["Data"]]$Value[which( ribogreenHigh732418[["Data"]]$set == "New" )] )

png( "C:/Data/BiomekPrep/densityHigh_732418.png" )
plotRibogreenDensities( densityOldHigh732418 , densityNewHigh732418 , study = 732418 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotHigh_732418.png" )
plotRibogreenScatterplot(  ribogreenHigh732418 , "Old" , "New" , study = 732418 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotHigh_732418.png" )
plotRibogreenBoxplot( ribogreenHigh732418 , study = 732418 )
dev.off()

ribogreenLow732418 <- makeRibogreen( layouts732418$oldLayoutLow , layouts732418$newLayoutLow , layouts732418$oldDataLow , layouts732418$newDataLow )

densityOldLow732418 <- density( ribogreenLow732418[["Data"]]$Value[which( ribogreenLow732418[["Data"]]$set == "Old" )] )
densityNewLow732418 <- density( ribogreenLow732418[["Data"]]$Value[which( ribogreenLow732418[["Data"]]$set == "New" )] )

png( "C:/Data/BiomekPrep/densityLow_732418.png" )
plotRibogreenDensities( densityOldLow732418 , densityNewLow732418 , study = 732418 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotLow_732418.png" )
plotRibogreenScatterplot(  ribogreenLow732418 , "Old" , "New" , study = 732418 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotLow_732418.png" )
plotRibogreenBoxplot( ribogreenLow732418 , study = 732418 )
dev.off()

# Study #732417

layouts732417 <- makeLayouts( layout , Data2 )

ribogreenHigh732417 <- makeRibogreen( layouts732417$oldLayoutHigh , layouts732417$newLayoutHigh , layouts732417$oldDataHigh , layouts732417$newDataHigh )

densityOldHigh732417 <- density( ribogreenHigh732417[["Data"]]$Value[which( ribogreenHigh732417[["Data"]]$set == "Old" )] , na.rm = TRUE )
densityNewHigh732417 <- density( ribogreenHigh732417[["Data"]]$Value[which( ribogreenHigh732417[["Data"]]$set == "New" )] , na.rm = TRUE )

png( "C:/Data/BiomekPrep/densityHigh_732417.png" )
plotRibogreenDensities( densityOldHigh732417 , densityNewHigh732417 , study = 732417 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotHigh_732417.png" )
plotRibogreenScatterplot(  ribogreenHigh732417 , "Old" , "New" , study = 732417 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotHigh_732417.png" )
plotRibogreenBoxplot( ribogreenHigh732417 , study = 732417 )
dev.off()

ribogreenLow732417 <- makeRibogreen( layouts732417$oldLayoutLow , layouts732417$newLayoutLow , layouts732417$oldDataLow , layouts732417$newDataLow )

densityOldLow732417 <- density( ribogreenLow732417[["Data"]]$Value[which( ribogreenLow732417[["Data"]]$set == "Old" )] , na.rm = TRUE )
densityNewLow732417 <- density( ribogreenLow732417[["Data"]]$Value[which( ribogreenLow732417[["Data"]]$set == "New" )] , na.rm = TRUE )

png( "C:/Data/BiomekPrep/densityLow_732417.png" )
plotRibogreenDensities( densityOldLow732417 , densityNewLow732417 , study = 732417 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotLow_732417.png" )
plotRibogreenScatterplot(  ribogreenLow732417 , "Old" , "New" , study = 732417 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotLow_732417.png" )
plotRibogreenBoxplot( ribogreenLow732417 , study = 732417 )
dev.off()

# Study #732537

layouts732537 <- makeLayouts( layout , Data3 )

ribogreenHigh732537 <- makeRibogreen( layouts732537$oldLayoutHigh , layouts732537$newLayoutHigh , layouts732537$oldDataHigh , layouts732537$newDataHigh )

densityOldHigh732537 <- density( ribogreenHigh732537[["Data"]]$Value[which( ribogreenHigh732537[["Data"]]$set == "Old" )] )
densityNewHigh732537 <- density( ribogreenHigh732537[["Data"]]$Value[which( ribogreenHigh732537[["Data"]]$set == "New" )] )

png( "C:/Data/BiomekPrep/densityHigh_732537.png" )
plotRibogreenDensities( densityOldHigh732537 , densityNewHigh732537 , study = 732537 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotHigh_732537.png" )
plotRibogreenScatterplot(  ribogreenHigh732537 , "Old" , "New" , study = 732537 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotHigh_732537.png" )
plotRibogreenBoxplot( ribogreenHigh732537 , study = 732537 )
dev.off()

ribogreenLow732537 <- makeRibogreen( layouts732537$oldLayoutLow , layouts732537$newLayoutLow , layouts732537$oldDataLow , layouts732537$newDataLow )

densityOldLow732537 <- density( ribogreenLow732537[["Data"]]$Value[which( ribogreenLow732537[["Data"]]$set == "Old" )] )
densityNewLow732537 <- density( ribogreenLow732537[["Data"]]$Value[which( ribogreenLow732537[["Data"]]$set == "New" )] )

png( "C:/Data/BiomekPrep/densityLow_732537.png" )
plotRibogreenDensities( densityOldLow732537 , densityNewLow732537 , study = 732537 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotLow_732537.png" )
plotRibogreenScatterplot(  ribogreenLow732537 , "Old" , "New" , study = 732537 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotLow_732537.png" )
plotRibogreenBoxplot( ribogreenLow732537 , study = 732537 )
dev.off()

# Study #732538

layouts732538 <- makeLayouts( layout , Data4 )

ribogreenHigh732538 <- makeRibogreen( layouts732538$oldLayoutHigh , layouts732538$newLayoutHigh , layouts732538$oldDataHigh , layouts732538$newDataHigh )

densityOldHigh732538 <- density( ribogreenHigh732538[["Data"]]$Value[which( ribogreenHigh732538[["Data"]]$set == "Old" )] )
densityNewHigh732538 <- density( ribogreenHigh732538[["Data"]]$Value[which( ribogreenHigh732538[["Data"]]$set == "New" )] )

png( "C:/Data/BiomekPrep/densityHigh_732538.png" )
plotRibogreenDensities( densityOldHigh732538 , densityNewHigh732538 , study = 732538 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotHigh_732538.png" )
plotRibogreenScatterplot(  ribogreenHigh732538 , "Old" , "New" , study = 732538 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotHigh_732538.png" )
plotRibogreenBoxplot( ribogreenHigh732538 , study = 732538 )
dev.off()

ribogreenLow732538 <- makeRibogreen( layouts732538$oldLayoutLow , layouts732538$newLayoutLow , layouts732538$oldDataLow , layouts732538$newDataLow )

densityOldLow732538 <- density( ribogreenLow732538[["Data"]]$Value[which( ribogreenLow732538[["Data"]]$set == "Old" )] )
densityNewLow732538 <- density( ribogreenLow732538[["Data"]]$Value[which( ribogreenLow732538[["Data"]]$set == "New" )] )

png( "C:/Data/BiomekPrep/densityLow_732538.png" )
plotRibogreenDensities( densityOldLow732538 , densityNewLow732538 , study = 732538 )
dev.off()

png( "C:/Data/BiomekPrep/scatterplotLow_732538.png" )
plotRibogreenScatterplot(  ribogreenLow732538 , "Old" , "New" , study = 732538 )
dev.off()

png( "C:/Data/BiomekPrep/boxplotLow_732538.png" )
plotRibogreenBoxplot( ribogreenLow732538 , study = 732538 )
dev.off()
