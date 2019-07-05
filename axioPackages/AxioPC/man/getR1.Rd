\name{getR1}
\alias{getR1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getR1(b0, b2, b1, b3, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b0}{ ~~Describe \code{b0} here~~ }
  \item{b2}{ ~~Describe \code{b2} here~~ }
  \item{b1}{ ~~Describe \code{b1} here~~ }
  \item{b3}{ ~~Describe \code{b3} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
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
function( b0 , b2 , b1 , b3 , p )
{
  p0 <- ( 1 - p )^2
  p1 <- 2 * p * (1 - p )
  p2 <- p^2
  logit <- c( b0 + b1 , b0 + b1 + b2 + b3 , b0 + b1 + 2 * ( b2 + b3 ) ) 
  prob <- exp( logit ) / ( 1 + exp( logit ) )
  return( prob \%*\% c( p0 , p1 , p2 ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
