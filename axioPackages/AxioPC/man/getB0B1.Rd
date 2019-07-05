\name{getB0B1}
\alias{getB0B1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getB0B1(b2, b3, p, r0, r1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b2}{ ~~Describe \code{b2} here~~ }
  \item{b3}{ ~~Describe \code{b3} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{r0}{ ~~Describe \code{r0} here~~ }
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
function( b2 , b3 , p , r0 , r1 )
{
  b0 <- bisecB0( a = -1000 , b = 1000 , b2 , epslon = 0.0001 , p = p , r0 = r0 )
  b1 <- bisecB1( a = -1000 , b = 1000 , b0 , b2 , b3 , epslon = 0.0001 , p = p , r1 = r1 )
  return( c( b0 , b1 ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
