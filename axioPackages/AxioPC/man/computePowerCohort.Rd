\name{computePowerCohort}
\alias{computePowerCohort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computePowerCohort(pg, pd, OR, N, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pg}{ ~~Describe \code{pg} here~~ }
  \item{pd}{ ~~Describe \code{pd} here~~ }
  \item{OR}{ ~~Describe \code{OR} here~~ }
  \item{N}{ ~~Describe \code{N} here~~ }
  \item{alpha}{ ~~Describe \code{alpha} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function( pg , pd , OR , N , alpha )
{
  p1 <- ( ( 1 + ( OR - 1 ) * ( pg + pd ) ) - sqrt( ( 1 + ( OR - 1 ) * ( pg + pd ) )^2 - 4 * pd * pg * ( OR - 1 ) * OR ) ) / ( 2 * pg * ( OR - 1 ) )
  p0 <- p1 / ( ( 1 - p1 ) * OR + p1 )
  ES <- sqrt( ( ( p1 - pd ) * pg )^2 / ( pd * pg ) + ( ( p0 - pd ) * ( 1 - pg ) )^2 / ( pd * ( 1 - pg ) ) + ( ( pd - p1 ) * pg )^2 / ( ( 1 - pd ) * pg ) + ( ( pd - p0 ) * ( 1 - pg ) )^2 / ( ( 1 - pd ) * ( 1 - pg ) ) )
  Pow.or <- OR
  for ( n in N )
  {
    pow <- pwr.chisq.test( w = ES , N = n , df = 1 , sigLevel = alpha )$power
    PowOr <- cbind( Pow.or , pow )
  }
  PowOr <- as.data.frame( PowOr )
  names( PowOr) <- c( "OR" , paste( "pow" , N , sep="." ) )
  return( PowOr )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
