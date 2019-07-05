\name{getNuTau2}
\alias{getNuTau2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getNuTau2}
\description{getNuTau2}
\usage{
getNuTau2(tau2, nu, sc, muNu, sigNu)
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
    w = 1
    m = 100
    z = getPostNuTau2(tau2, nu, sc, muNu, sigNu) - rexp(1)
    u = runif(1)
    L = nu - w * u
    R = L + w
    v = runif(1)
    J = floor(m * v)
    K = (m - 1) - J
    L = max(0, L)
    while (J > 0 && L > 0 && z < getPostNuTau2(tau2, L, sc, muNu, 
        sigNu)) {
        L = L - w
        L = max(0, L)
        J = J - 1
    }
    while (K > 0 && z < getPostNuTau2(tau2, nu, R, muNu, sigNu)) {
        R = R + w
        K = K - 1
    }
    u = runif(1)
    newParam = L + u * (R - L)
    while (z > getPostNuTau2(tau2, newParam, sc, muNu, sigNu)) {
        if (newParam < nu) {
            L = newParam
        }
        else {
            R = newParam
        }
        u = runif(1)
        newParam = L + u * (R - L)
    }
    return(newParam)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
