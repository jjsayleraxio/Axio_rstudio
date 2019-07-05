\name{getPostNuTau2}
\alias{getPostNuTau2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getPostNuTau2}
\description{getPostNuTau2}
\usage{
getPostNuTau2(tau2, nu, sc, muNu, sigNu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau2}{
%%     ~~Describe \code{tau2} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{sc}{
%%     ~~Describe \code{sc} here~~
}
  \item{muNu}{
%%     ~~Describe \code{muNu} here~~
}
  \item{sigNu}{
%%     ~~Describe \code{sigNu} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tau2, nu, sc, muNu, sigNu) 
{
    logPrior = dgamma(nu, muNu, sigNu, log = TRUE)
    logLike = sum(get.dichi2.log(tau2, nu, sc))
    return(logPrior + logLike)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
