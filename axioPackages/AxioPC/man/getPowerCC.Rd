\name{getPowerCC}
\alias{getPowerCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPowerCC(ORg, ORgxt, N00, N01, N10, N11, p, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ORg}{ ~~Describe \code{ORg} here~~ }
  \item{ORgxt}{ ~~Describe \code{ORgxt} here~~ }
  \item{N00}{ ~~Describe \code{N00} here~~ }
  \item{N01}{ ~~Describe \code{N01} here~~ }
  \item{N10}{ ~~Describe \code{N10} here~~ }
  \item{N11}{ ~~Describe \code{N11} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
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
function( ORg , ORgxt , N00 , N01 , N10 , N11 , p , alpha )
{
  p0 <- ( 1 - p )^2
  p1 <- 2 * p * ( 1 - p )
  p2 <- p^2
  b2 <- log( ORg )
  b3 <- log( ORgxt )
  N0 <- N00 + N01
  N1 <- N11 + N10
  r0 <- N01 / N0
  r1 <- N11 / N1
  N <- N0 + N1
  b01 <- getB0B1( b2 , b3 , p , r0 , r1 )
  b0 <- b01[1]
  b1 <- b01[2]
  logit <- c( b0 , b0 + b2 , b0 + 2 * b2 , b0 + b1 , b0 + b1 + b2 + b3 , b0 + b1 + 2 * ( b2 + b3 ) ) 
  temp <- exp( -logit ) / ( ( 1 + exp( -logit ) )^2 )
  n.vec <- c( c( p0 , p1 , p2 ) * N0 , c( p0 , p1 , p2 ) * N1 ) 

  fisher <- matrix( NA , ncol = 4 , nrow = 4 )
  fisher[1,1] <- temp \%*\% n.vec
  fisher[1,2] <- fisher[2,1] <- temp[4:6] \%*\% n.vec[4:6]
  fisher[1,3] <- fisher[3,1] <- temp[2] * n.vec[2] + temp[5] * n.vec[5] + 2 * 
                 ( temp[3] * n.vec[3] + temp[6] * n.vec[6] )
  fisher[1,4] <- fisher[4,1] <- temp[5] * n.vec[5] + 2 * temp[6] * n.vec[6]
  fisher[2,2] <- temp[4:6] \%*\% n.vec[4:6]
  fisher[2,3] <- fisher[3,2] <- temp[5] * n.vec[5] + 2 * temp[6] * n.vec[6]
  fisher[2,4] <- fisher[4,2] <- temp[5] * n.vec[5] + 2 * temp[6] * n.vec[6]
  fisher[3,3] <- temp[2] * n.vec[2] + temp[5] * n.vec[5] + 4 * ( temp[3] * n.vec[3] + temp[6] * n.vec[6] )
  fisher[3,4] <- fisher[4,3] <- temp[5] * n.vec[5] + 4 * temp[6] * n.vec[6]
  fisher[4,4] <- temp[5] * n.vec[5] + 4 * temp[6] * n.vec[6]

  se <- sqrt( solve( fisher )[ 4 , 4 ] )
  cutoff <- abs( qnorm( alpha ) ) - b3 / se
  return( 1 - pnorm( cutoff ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
