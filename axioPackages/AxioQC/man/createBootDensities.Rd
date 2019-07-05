\name{createBootDensities}
\alias{createBootDensities}
\title{
Create bootstrap densities.
}
\description{
Function to create the bootstrap densities for plotting and finding the limits.
}
\usage{
createBootDensities( data , date = "Date" , level = "100ng" , type = "quarterly" , nsamples = 20000 , pCutoff = 0.01 , subset = NULL , bootD = TRUE , cumulative = FALSE )
}
\arguments{
  \item{data}{ list, data object.  A \code{list} object with the bootstrap density data contained inside. }
  \item{date}{ character string, name of date column in data. }
  \item{level}{ character string, column name of data for density plot. }
  \item{type}{ character string, type of date segmentation.  Must be one of: "quarterly", "bimonthly", "monthly", "biweekly", "weekly", or "daily". }
  \item{nsamples}{ numeric, Number of bootstrap samples to generate. }
  \item{pCutoff}{ One tail limit for selecting outer boundary for density limits. }
  \item{subset}{ vector, index of subset of data if requested. }
  \item{bootD}{ logical, keep the bootstrap densities. }
  \item{cumulative}{ logical, are the densities cumulative. }
}
\details{
Creates a list of bootstrap densities.
}
\references{Kernel Smoothing to Improve Bootstrap Confidence Intervals (J. Roy. Statist. Soc., B, 59 (1997), 821-838)}
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
