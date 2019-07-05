\name{getSampleElstonPGxNormalVec}
\alias{getSampleElstonPGxNormalVec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSampleElstonPGxNormalVec(p = 0.2, gmeans, sigma2, r = 0.5, gaction = c("additive", "dominant", "recessive"), alpha = 0.05, beta = 0.8, type = c("main", "interaction"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{gmeans}{ ~~Describe \code{gmeans} here~~ }
  \item{sigma2}{ ~~Describe \code{sigma2} here~~ }
  \item{r}{ ~~Describe \code{r} here~~ }
  \item{gaction}{ ~~Describe \code{gaction} here~~ }
  \item{alpha}{ ~~Describe \code{alpha} here~~ }
  \item{beta}{ ~~Describe \code{beta} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
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
function (p = 0.2, gmeans, sigma2, r = 0.5, gaction = c("additive", 
    "dominant", "recessive"), alpha = 0.05, beta = 0.8, type = c("main", 
    "interaction")) 
{
    return(as.vector(lapply(gmeans, getSampleElstonPGxNormal, 
        sigma2, p, r, gaction, alpha, beta, type)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
