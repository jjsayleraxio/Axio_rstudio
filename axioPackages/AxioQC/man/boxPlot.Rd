\name{boxPlot}
\alias{boxPlot}
\title{
Function to plot QC summaries
}
\description{
Plots boxplots by segmentation and time.
}
\usage{boxPlot( data , x = "x" , y = "y" , time , xlab , ylab , main , timeName )
}
\arguments{
  \item{data}{ data object. }
  \item{x}{ character string, name of x axis variable. }
  \item{y}{ character string, type of y axis variable. }
  \item{time}{ variable, column name of data for time. }
  \item{xlab}{ character string, label of x axis. }
  \item{ylab}{ character string, label of y axis. }
  \item{main}{ character string, title of plot. }
  \item{timeName}{ character string, label of time legend. }
}
\details{
Creates a box plot by time and segmentation.
}
\examples{
data( x10NoOutliers )

boxPlot( data = x10NoOutliers , x = "Operator" , y = "Concentration" , time = "TimeOfDay" , xlab = "Operator" , ylab = "Concentration" , main = "Boxplot of Concentration by Operator and Time of Day" , timeName = "Time of Day" )
boxPlot( data = x10NoOutliers , x = "RobotID" , y = "Concentration" , time = "TimeOfDay" , xlab = "Robot ID" , ylab = "Concentration" , main = "Boxplot of Concentration by Operator and Time of Day" , timeName = "Time of Day" )
boxPlot( data = x10NoOutliers , x = "MultiRobotFlag" , y = "Concentration" , time = "TimeOfDay" , xlab = "Multiple Robot Flag" , ylab = "Concentration" , main = "Boxplot of Concentration by Operator and Time of Day" , timeName = "Time of Day" )
boxPlot( data = x10NoOutliers , x = "MultiOperatorFlag" , y = "Concentration" , time = "TimeOfDay" , xlab = "Multiple Operator Flag" , ylab = "Concentration" , main = "Boxplot of Concentration by Operator and Time of Day" , timeName = "Time of Day" )
boxPlot( data = x10NoOutliers , x = "RedoFlag" , y = "Concentration" , time = "TimeOfDay" , xlab = "Redo Flag" , ylab = "Concentration" , main = "Boxplot of Concentration by Operator and Time of Day" , timeName = "Time of Day" )
}
\keyword{boxplot}
