\name{getCountData}
\alias{getCountData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Count Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCountData(x, multiplexed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{multiplexed}{
%%     ~~Describe \code{multiplexed} here~~
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
##--  or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, multiplexed = FALSE) 
{
    if (multiplexed) {
        xNames <- unlist(strsplit(names(x$sampleAttributes)[2], 
            " ", fixed = TRUE))
    }
    else {
        xNames <- names(x$sampleAttributes)[2]
    }
    x <- x$codeSummary
    cMat <- matrix(nrow = length(unique(x$Name)), ncol = length(xNames))
    rownames(cMat) <- unique(x$Name)
    colnames(cMat) <- xNames
    controlIndex <- which(x[, "CodeClass"] \%in\% c("Positive", 
        "Negative"))
    controls <- x[controlIndex, "Count"]
    names(controls) <- x[controlIndex, "Name"]
    x <- x[-controlIndex, ]
    x$Group <- as.integer(substr(x$CodeClass, nchar(as.character(x$CodeClass)), 
        nchar(as.character(x$CodeClass))))
    for (i in seq(1, length(xNames))) {
        indx <- which(x$Group \%in\% i)
        cMat[c(names(controls), as.character(x[indx, "Name"])), 
            i] <- c(controls, x[indx, "Count"])
    }
    return(cMat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
