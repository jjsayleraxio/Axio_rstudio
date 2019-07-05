\name{plot.powerObject}
\alias{plot.powerObject}
\title{Plot Power Calculation Values}
\description{
Generates a series of plots of power calculation statistics.  The upper plot is power versus minor allele frequency (MAF), the middle plot is sample size versus MAF, and the lower set of plots are the genetic effects tested.
}
\usage{
plot.powerObject( x , y = NULL , SD = NULL , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The power object}
  \item{y}{Unused at present}
  \item{SD}{Alternative description of standard deviation for when each effect has different variability.}
  \item{\dots}{Unused at present}
}
\details{
Create a set of plots of power calculation values.
}
\value{
Returns a ggplot2 object of the set of plots.
}
\references{None}
\author{David A. Henderson <\email{David.Henderson@Covance.Com}>}
\note{
Still a work in progress.  API may change without notice.
}
\seealso{\code{\link{computePowerElstonPGx}},\code{\link{getPowerElstonPGxNormal}},\code{\link{getSampleElstonPGxNormal}}}
\examples{
#Main Effects
geneticMeans <- list( additive = c( 100 , 110 , 120 ) , dominance = c( 115 , 115 , 120 ) )
geneticVariance <- 3^2
powerElston <- computePowerElstonPGx( geneticMeans , geneticVariance , p = seq( 0.05 , 0.2 , 0.05 ) , n = 100 , r = 0.5 ,
                          gaction = "additive" , alpha = 0.05 , beta = 0.8 , type = "main" )
plot( powerElston )
# Interaction
geneticMeans <- list( additive = matrix( c( 100 , 112 , 118 , 120 , 117 , 108 ) , nrow = 2 , byrow = TRUE ) ,
                      dominance = matrix( c( 100 , 117 , 115 , 117 , 115 , 105 ) , nrow = 2 , byrow = TRUE ) )
powerElston <- computePowerElstonPGx( geneticMeans , geneticVariance , p = seq( 0.05 , 0.2 , 0.05 ) , n = 100 ,
                      r = c( 0.5 , 0.4 ) , gaction = "additive" , alpha = 0.05 , beta = 0.8 , type = "interaction" )
plot( powerElston )
}
\keyword{power calculation}
\keyword{sample size}
