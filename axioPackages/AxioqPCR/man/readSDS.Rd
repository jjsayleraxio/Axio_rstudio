\name{readSDS}
\alias{readSDS}

\title{
Function to read in qPCR SDS files
}
\description{
Reads in SDS files for use in combining multiple files in to one.
}
\usage{
readSDS(dataFile = "", qPCRPlatform = "STD", plate = NULL)
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
sdsData <- readSDS( dataFile = "sds.txt" , qPCRPlatform = "STD" , plate = 1 ) 
}
\keyword{ qPCR }
