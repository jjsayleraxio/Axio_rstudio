library( XLConnect )

source( "c:/Documents and Settings/dhender1/workspace/qPCR/R/DoBy.R" )

dataWB <- loadWorkbook( "C:/Data/BiomekPrep/final_concentration_comparisons_eaw_4Oct2011.xls" )

Data <- readWorksheet( dataWB , sheet = "data" , startRow = 1 , endRow = 726 , startCol = 1 , endCol = 13 , header = TRUE )

names( Data ) <- c( "PlateBC" , "WellIndex1" , "Barcode" , "Volume1" , "Concentration1" , "Concentration2" , "Volume2" , "WellIndex2" , "Diff" , "Avg" , "PctDiff" , "abs(PctDiff)" , "Row" )

clrs <- rainbow( 94 )

clrmap <- vector( length = nrow( Data ) )
i <- 1
for ( j in unique( Data$WellIndex1 ) )
{
  clrmap[which( Data$WellIndex1 == j )] <- clrs[i]
  i <- i + 1
}

plot( Data$Concentration1 , Data$Concentration2 , col = clrmap )
plot( Data$Volume1 , Data$Volume2 , col = clrmap )

clrs <- terrain.colors( 8 )

clrmap <- vector( length = nrow( Data ) )
i <- 1
for ( j in unique( Data$Row ) )
{
  clrmap[which( Data$Row == j )] <- clrs[i]
  i <- i + 1
}

pdf( "c:/Documents and Settings/dhender1/My Documents/Biomek/final_concentration.pdf" )
plot( Data$Concentration1 , Data$Concentration2 , col = clrmap , xlab = "Concentration 1" , ylab = "Concentration 2" , main = "Concentration 1 by Concentration 2" )
plot( Data$Volume1 , Data$Volume2 , col = clrmap , xlab = "Volume 1" , ylab = "Volume 2" , main = "Volume 1 by Volume 2" )
boxplot( Concentration1 ~ Row , Data , main = "Boxplot of Concentration 1 by Plate Row" )
boxplot( Concentration2 ~ Row , Data , main = "Boxplot of Concentration 2 by Plate Row" )
boxplot( Volume1 ~ Row , Data , main = "Boxplot of Volume 1 by Plate Row" )
boxplot( Volume2 ~ Row , Data , main = "Boxplot of Volume 2 by Plate Row" )

summaryData1 <- DoBy( "Concentration1" , c( "Row" ) , data = Data , FUN = "sd" , na.rm = TRUE )
summaryData1$color <- unlist( apply( summaryData1 , 1 , function( x ) { ifelse( x[2] < 30 , "blue" , "red" ) } ) )
summaryData2 <- DoBy( "Concentration2" , c( "Row" ) , data = Data , FUN = "sd" , na.rm = TRUE )
summaryData2$color <- unlist( apply( summaryData2 , 1 , function( x ) { ifelse( x[2] < 30 , "green" , "orange" ) } ) )
ylimit <- c( min( c( summaryData1$Concentration1.sd , summaryData2$Concentration2.sd ) ) , max( c( summaryData1$Concentration1.sd , summaryData2$Concentration2.sd ) ) )
ylimit <- ylimit + c( -0.05 , 0.10 ) * ylimit
plot( summaryData1$Row , summaryData1$Concentration1.sd , xlab = "Row" , ylab = "Standard Deviation" , col = summaryData1$color , main = "Variability by Row" , ylim = ylimit )
points( summaryData2$Row , summaryData2$Concentration2.sd , pch = 2 , col = summaryData2$color )
legend( "topleft" , legend = c( "Run 1, Low Variability" , "Run 1, High Variability" ) , col = c( "blue" , "red" ) , pch = 1 )
legend( "topright" , legend = c( "Run 2, Low Variability" , "Run 2, High Variability" ) , col = c( "green" , "orange" ) , pch = 2 )
dev.off()
