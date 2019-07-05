\name{print.powerObject}
\alias{print.powerObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.powerObject(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
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
function (object) 
{
    print("Effect sizes examined")
    print(object@effects)
    print("Minor allele frequencies examined")
    print(object@MAF)
    print("Standard Deviation examined")
    print(object@SD)
    print("N examined")
    print(object@N)
    print("Type I Error examined")
    print(object@alpha)
    print("Power examined")
    print(object@beta)
    print("Proportion in treatment arm")
    print(object@proportion)
    if (!is.null(object@power)) {
        print("Power from calculations")
        tmp <- as.data.frame(object@power)
        rownames(tmp) <- object@effects
        colnames(tmp) <- object@MAF
        print(tmp)
    }
    if (!is.null(object@sampleSize)) {
        print("Sample Size from calculations")
        tmp <- as.data.frame(object@sampleSize)
        rownames(tmp) <- object@effects
        colnames(tmp) <- object@MAF
        print(tmp)
    }
    invisible(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
