\name{smoothPlot}
\alias{smoothPlot}
\title{
Function to plot QC summaries
}
\description{
Plots boxplots and smoothed lines by segmentation and time.
}
\usage{
smoothPlot( data , x = "x" , y = "y" , wQCDate , time , xlab , ylab , xBreakMajor , xBreakMinor , main , timeName )
}
\arguments{
  \item{data}{ data object }
  \item{x}{ character string, name of x axis variable. }
  \item{y}{ character string, type of y axis variable. }
  \item{wQCDate}{ variable, time segmentation variable. }
  \item{time}{ variable, column name of data for time. }
  \item{xlab}{ character string, label of x axis. }
  \item{ylab}{ character string, label of y axis. }
  \item{xBreakMajor}{ character string, time segmentation for major breaks. }
  \item{xBreakMinor}{ character string, time segmentation for minor breaks. }
  \item{main}{ character string, title of plot. }
  \item{timeName}{ character string, label of time legend. }
}
\details{
Creates a box plot by time and segmentation.
}
\examples{
data( x10NoOutliers )

smoothPlot( x10NoOutliers , "ExecutionDate" , "Concentration" , "wExecutionDate" , "TimeOfDay" , "Execution Date" , "Concentration" , "2 weeks" , "1 week" , "Boxplot of Concentration by Execution Date and Time of Day" , "Time of Day" )
smoothPlot( x10NoOutliers , "QCDate" , "Concentration" , "wExecutionDate" , "TimeOfDay" , "QC Date" , "Concentration" , "2 weeks" , "1 week" , "Boxplot of Concentration by QC Date and Time of Day" , "Time of Day" )
smoothPlot( x10NoOutliers , "QCDate" , "Volume" , "wExecutionDate" , "TimeOfDay" , "QC Date" , "Volume" , "2 weeks" , "1 week" , "Boxplot of Volume by QC date and Time of Day" , "Time of Day" )
}
\keyword{boxplot}
