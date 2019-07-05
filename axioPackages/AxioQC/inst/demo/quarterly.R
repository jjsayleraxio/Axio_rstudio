library( XLConnect )

source( "c:/Documents and Settings/dhender1/workspace/BiomekPrep/smoothedBootstrap.R" )
source( "c:/Documents and Settings/dhender1/workspace/BiomekPrep/helperFunctions.R" )

dataWB <- loadWorkbook( "C:/Data/BiomekPrep/DNA_Controls_by_Quarter.xls" )

q42010 <- readWorksheet( dataWB , sheet = "Q4 2010" , header = TRUE )
q12011 <- readWorksheet( dataWB , sheet = "Q1 2011" , header = TRUE )
q22011 <- readWorksheet( dataWB , sheet = "Q2 2011" , header = TRUE )
q32011 <- readWorksheet( dataWB , sheet = "Q3 2011" , header = TRUE )

q42010[,1] <- as.Date( q42010[,1] )
q12011[,1] <- as.Date( q12011[,1] )
q22011[,1] <- as.Date( q22011[,1] )

q32011 <- q32011[-nrow( q32011 ),]
q32011[15,2] <- "Fri Aug 12 00:00:00 PDT 2011"
q32011[,2] <- as.Date( paste( substr( q32011[,2] , 9 , 10 ) , substr( q32011[,2] , 5 , 7 ) , substr( q32011[,2] , 25 , 28 ) , sep = "-" ) , format = "%d-%B-%Y" )

names( q42010 ) <- c( "Date" , "100ng" , "10ng" )
names( q12011 ) <- c( "Date" , "100ng" , "10ng" )
names( q22011 ) <- c( "Date" , "100ng" , "10ng" )
names( q32011 ) <- c( "Run" , "Date" , "100ng" , "10ng" )

alldata <- rbind( q42010 , q12011 , q22011 , q32011[,-1] )
if ( any( alldata[,"100ng"] == 0 , na.rm = TRUE ) )
{
  alldata[which( alldata[,"100ng"] == 0 ),] <- NA
}
if ( any( alldata[,"10ng"] == 0 , na.rm = TRUE ) )
{
  alldata[which( alldata[,"10ng"] == 0 ),] <- NA
}
if ( any( alldata[,"100ng"] > 150 | alldata[,"100ng"] < 50 , na.rm = TRUE ) )
{
  alldata[which( alldata[,"100ng"] > 150 | alldata[,"100ng"] < 50 ),] <- NA
}
if ( any( alldata[,"10ng"] > 40 | alldata[,"10ng"] < 0 , na.rm = TRUE ) )
{
  alldata[which( alldata[,"10ng"] > 40 | alldata[,"10ng"] < 0 ),] <- NA
}

# Sequential stuff

allDataList10 <- createBootDensities( alldata , level = "10ng" )
allDataList100 <- createBootDensities( alldata , level = "100ng" )

plotDensities( allDataList100 , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , PNG = TRUE )
plotDensities( allDataList100 , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q4-2010" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList100 , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q1-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList100 , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q2-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList100 , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q3-2011" , plotBootD = TRUE ,
    PNG = TRUE )

plotDensities( allDataList10 , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , PNG = TRUE )
plotDensities( allDataList10 , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q4-2010" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList10 , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q1-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList10 , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q2-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList10 , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , subset = "Q3-2011" , plotBootD = TRUE ,
    PNG = TRUE )

plotLimitWidths( allDataList10 , pathRoot = "c:/Data/BiomekPrep/" , level = "10ng" , PNG = TRUE )
plotLimitWidths( allDataList100 , pathRoot = "c:/Data/BiomekPrep/" , level = "100ng" , PNG = TRUE )
    
# Cumulative stuff

allDataList10c <- createBootDensities( alldata , level = "10ng" , cumulative = TRUE )
allDataList100c <- createBootDensities( alldata , level = "100ng" , cumulative = TRUE )

plotDensities( allDataList100c , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , PNG = TRUE )
plotDensities( allDataList100c , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q4-2010" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList100c , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q1-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList100c , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q2-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList100c , level = "100ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q3-2011" , plotBootD = TRUE ,
    PNG = TRUE )

plotDensities( allDataList10c , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , PNG = TRUE )
plotDensities( allDataList10c , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q4-2010" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList10c , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q1-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList10c , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q2-2011" , plotBootD = TRUE ,
    PNG = TRUE )
plotDensities( allDataList10c , level = "10ng" , pathRoot = "c:/Data/BiomekPrep/" ,
    main = "Density by Quarter" , cumulative = TRUE , subset = "Q3-2011" , plotBootD = TRUE ,
    PNG = TRUE )

plotLimitWidths( allDataList10c , pathRoot = "c:/Data/BiomekPrep/" , level = "10ng" , cumulative = TRUE , PNG = TRUE )
plotLimitWidths( allDataList100c , pathRoot = "c:/Data/BiomekPrep/" , level = "100ng" , cumulative = TRUE , PNG = TRUE )
