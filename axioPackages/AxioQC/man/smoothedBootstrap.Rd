\name{smoothedBootstrap}
\alias{smoothedBootstrap}
\title{
Perform smoothed bootstrap
}
\description{
Performs a smoothed bootstrap for use in small data sets
}
\usage{
smoothedBootstrap( x , n )
}
\arguments{
  \item{x}{ vector, data object. }
  \item{n}{ integer, number of bootstrap samples. }
}
\details{
Performs a smoothed bootstrap.
}
\references{Kernel Smoothing to Improve Bootstrap Confidence Intervals (J. Roy. Statist. Soc., B, 59 (1997), 821-838)}
\examples{
data( allData )

allDataBoot <- smoothedBootstrap( allData[which( !is.na( allData[,"10ngul"] ) ),"10ngul"] , n = 1000 )
}
\keyword{bootstrap}
\keyword{smoothing}
