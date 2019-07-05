\name{plotLimitWidths}
\alias{plotLimitWidths}
\title{
Plot the limit widths and a slope around them over time
}
\description{
Plot the limit widths and a slope around them over time
}
\usage{
plotLimitWidths( data , pathRoot = "./" , level = "100ng" , type = "quarterly" , main = "Interval Width by Quarter" , xlab = "Quarter" , ylab = "Interval Width" , cumulative = FALSE , PNG = FALSE )
}
\arguments{
  \item{data}{ list, data object.  A \code{list} object with the bootstrap density data contained inside. }
  \item{pathRoot}{ character string, root folder to place png plot in if requested. }
  \item{level}{ character string, column name of data for density plot. }
  \item{type}{ character string, type of date segmentation.  Must be one of: "quarterly", "bimonthly", "monthly", "biweekly", "weekly", or "daily". }
  \item{main}{ character string, title of plot. }
  \item{xlab}{ character string, label of x axis. }
  \item{ylab}{ character string, label of y axis. }
  \item{cumulative}{ logical, are the densities cumulative. }
  \item{PNG}{ logical, create png file. }
}
\details{
Plot the limit widths and a slope around them over time.
}
\examples{
data( allData )

allDataList10 <- createBootDensities( allData , level = "10ngul" )
allDataList10w <- createBootDensities( allData , level = "10ngul" , type = "monthly" )
allDataList100 <- createBootDensities( allData , level = "100ngul" )

( plotLimitWidths( allDataList10 , level = "10ngul" , main = "Limit Widths by Quarter" ) )
( plotLimitWidths( allDataList10w , level = "10ngul" , type = "monthly" , main = "Limit Widths by Month" ) )
( plotLimitWidths( allDataList100 , level = "100ngul" , main = "Limit Widths by Quarter" ) )
}
\keyword{limits}
