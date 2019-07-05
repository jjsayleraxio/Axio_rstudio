\name{plotDensities}
\alias{plotDensities}
\title{ plotDensities: Function to plot bootstrap densities }
\description{
Plots line densities by date segment of bootstrap samples.
}
\usage{plotDensities( data , date = "Date" , type = "quarterly" , level = "100ng" , pathRoot = "./" , xlab = "Concentration" , ylab = "Density" , main = "Density by Quarter, 100ng" , subset = NULL , showLimits = NULL , plotBootD = FALSE , cumulative = FALSE , PNG = FALSE )
}
\arguments{
  \item{data}{ list, data object.  A \code{list} object with the bootstrap density data contained inside. }
  \item{date}{ character string, name of date column in data. }
  \item{type}{ character string, type of date segmentation.  Must be one of: "quarterly", "bimonthly", "monthly", "biweekly", "weekly", or "daily". }
  \item{level}{ character string, column name of data for density plot. }
  \item{pathRoot}{ character string, root folder to place png plot in if requested. }
  \item{xlab}{ character string, label of x axis. }
  \item{ylab}{ character string, label of y axis. }
  \item{main}{ character string, title of plot. }
  \item{subset}{ vector, index of subset of data if requested. }
  \item{showLimits}{ vector of segments to display. }
  \item{plotBootD}{ logical, plot the bootstrap densities. }
  \item{cumulative}{ logical, are the densities cumulative. }
  \item{PNG}{ logical, create png file. }
}
\details{
Creates a density plot with lines for each calendar segmentation listed in the data.
}
\examples{
data( allData )

allDataList10 <- createBootDensities( allData , level = "10ngul" )
allDataList10w <- createBootDensities( allData , level = "10ngul" , type = "monthly" )
allDataList100 <- createBootDensities( allData , level = "100ngul" )

( plotDensities( allDataList10 , level = "10ngul" , main = "Density by Quarter" , showLimits = c( "Q2-2012" , "Q3-2012" , "Q4-2012" , "Q1-2013" , "Q2-2013" , "All" ) ) )
( plotDensities( allDataList10w , level = "10ngul" , type = "monthly" , main = "Density by Month" , showLimits = c( "Q2-2012" , "Q3-2012" , "Q4-2012" , "Q1-2013" , "Q2-2013" , "All" ) ) )
( plotDensities( allDataList100 , level = "100ngul" , main = "Density by Quarter" , showLimits = c( "Q2-2012" , "Q3-2012" , "Q4-2012" , "Q1-2013" , "Q2-2013" , "All" ) ) )
}
\keyword{bootstrap}
\keyword{density}
