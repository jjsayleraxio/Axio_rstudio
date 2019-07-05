\name{getPostSigTau2}
\alias{getPostSigTau2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getPostSigTau2}
\description{getPostSigTau2}
\usage{
getPostSigTau2(tau2, nu, sc, s0, muSc, sigSc)
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
  \item{s0}{
%%     ~~Describe \code{s0} here~~
}
  \item{muSc}{
%%     ~~Describe \code{muSc} here~~
}
  \item{sigSc}{
%%     ~~Describe \code{sigSc} here~~
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
function (tau2, nu, sc, s0, muSc, sigSc) 
{
    logPrior = dgamma(sc, muSc, sigSc, log = TRUE)
    logLike = sum(get.dichi2.log(tau2, nu, sc + s0))
    return(logPrior + logLike)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
