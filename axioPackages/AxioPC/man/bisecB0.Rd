\name{bisecB0}
\alias{bisecB0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bisecB0(a, b, b2, epslon, p, r0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ ~~Describe \code{a} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{b2}{ ~~Describe \code{b2} here~~ }
  \item{epslon}{ ~~Describe \code{epslon} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{r0}{ ~~Describe \code{r0} here~~ }
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
function( a , b , b2 , epslon , p , r0 )
{ 
  while( abs( b - a ) > epslon )
  {
    const <- ( a + b ) / 2
    fa <- getR0( b0 = a , b2 = b2 , p ) - r0
    fc <- getR0( b0 = const , b2 = b2 , p ) - r0
    if ( fa * fc <= 0 )
    {
      b <- const
    }
    if ( fa * fc > 0 )
    {
      a <- const
    }
  }
  return( ( a + b ) / 2 )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
