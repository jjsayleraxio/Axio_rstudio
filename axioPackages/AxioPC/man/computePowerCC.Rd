\name{computePowerCC}
\alias{computePowerCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computePowerCC(pg, pd, OR, N, alpha)
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
  pgD <- p1 * pg / pd
  pgND <- ( 1 - p1 ) * pg / ( 1 - pd )
  ES <- 2 * asin( sqrt( pgD ) ) - 2 * asin( sqrt( pg.ND ) )
  PowOR <- OR
  if ( !is.null( dim( N ) ) )
  {
    for( k in 1:dim( N )[1] )
    {
      print( k )
      nCase <- N[k,1]
      nControl <- N[k,2]

      pow <- pwr.2p2n.test( ES , n1 = n.case , n2 = n.control , sig.level = alpha , alternative = "two.sided" )$power
      PowOR <- cbind( PowOR , pow )
    }
    PowOR <- as.data.frame( PowOR )
    names( PowOR) <- c( "OR ", paste( "pow" , N[,1] , N[,2] , sep="." ) )
  }
  else
  {
    pow <- pwr.2p2n.test( ES , n1 = N[1] , n2 = N[2] , sig.level = alpha , alternative = "two.sided" )$power
    PowOR <- cbind( PowOR , pow )
    PowOR <- as.data.frame( PowOR )
    names( PowOR ) <- c( "OR" , paste( "pow" , N[1] , N[2] , sep="." ) )
  }
  return( PowOR )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
