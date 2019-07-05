\name{getSample}
\alias{getSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSample(sy1, sy0, p, b2, b3, alpha, power, r1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sy1}{ ~~Describe \code{sy1} here~~ }
  \item{sy0}{ ~~Describe \code{sy0} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{b2}{ ~~Describe \code{b2} here~~ }
  \item{b3}{ ~~Describe \code{b3} here~~ }
  \item{alpha}{ ~~Describe \code{alpha} here~~ }
  \item{power}{ ~~Describe \code{power} here~~ }
  \item{r1}{ ~~Describe \code{r1} here~~ }
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
function( sy1 , sy0 , p , b2 , b3 , alpha , power , r1 )
{
  r0 <- 1 - r1
  var1 <- sy1^2
  var0 <- sy0^2
  sigma2 <- r1 * var1 + r0 * var0 - 2 * p * ( 1 - p ) * ( b2^2 + 2 * r1 * b2 * b3 + r1 * b3^2 )
  const <- abs( qnorm( 1 - alpha / 2 ) )
  b <- qnorm( power )
  return( ( b + const ) * ( b + const ) * sigma2 / ( r * ( 1 - r ) * b3^2 * 2 * p * ( 1 - p ) ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
