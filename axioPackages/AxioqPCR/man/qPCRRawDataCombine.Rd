\name{qPCRRawDataCombine}
\alias{qPCRRawDataCombine}
\title{
Function to combine multiple SDS files into one
}
\description{
Function that takes either a directory or a list of files and a directory and combines the SDS files therein into a single file.
}
\usage{
qPCRRawDataCombine(sdsDirectory = getwd(), outputFileName = "", qPCRPlatform = c("STD", "TLDA"), fileList = NULL, excel = FALSE, thresholdLimit = 0.2 , undetermined = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdsDirectory}{
Directory containing the data files
}
  \item{outputFileName}{
Name of the output file to be created
}
  \item{qPCRPlatform}{
Either STD or TLDA
}
  \item{fileList}{
Default is NULL which then inputs all files with a .txt extening in the sdsDirectory.  Otherwise a list of the file names in the sdsDirectory which are to be input.
}
  \item{excel}{
Create an Excel spreadsheet.  Default is FALSE.
}
  \item{thresholdLimit}{
Set threshold limit.  Default is 0.2.
}
  \item{undetermined}{
Numeric value to set undetermined values to.  Default is NULL which replaces undetermined values with NA.
}
}
\details{
Function that takes either a directory or a list of files and a directory and combines the SDS files therein into a single file.
}
\value{
Creates a file with the data concatenated into a single file.
}
\author{
David Henderson
}
\examples{
qPCRRawDataCombine( sdsDirectory = getwd() , outputFileName = "output.csv" , qPCRPlatform = "STD" , 
    fileList = NULL , excel = FALSE , thresholdLimit = 0.2 )
qPCRRawDataCombine( sdsDirectory = getwd() , outputFileName = "output.xlsx" , qPCRPlatform = "STD" , 
    fileList = NULL , excel = TRUE , thresholdLimit = 0.2 , undetermined = 38 )
}
\keyword{ qPCR }
