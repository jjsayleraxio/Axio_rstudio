\name{readSDStxt}
\alias{readSDStxt}

\title{
Function to read in qPCR SDS files
}
\description{
Reads in SDS files for use in combining multiple files in to one.
}
\usage{
readSDStxt(dataFile = "", qPCRPlatform = "STD", plate = NULL, undetermined = NULL)
}
\arguments{
  \item{dataFile}{
SDS data file name
}
  \item{qPCRPlatform}{
STD or TLDA
}
  \item{plate}{
Plate number
}
  \item{undetermined}{
Numeric value to set undetermined values to.  Default is NULL which replaces undetermined values with NA.
}
}
\details{
Reads in SDS files for use in combining multiple files in to one.
}
\value{
returns a data frame with the values from the input file.
}
\author{
David Henderson
}
\examples{
sdsData <- readSDStxt( dataFile = "sds.txt" , qPCRPlatform = "STD" , plate = 1, undetermined = 38 ) 
}
\keyword{ qPCR }
