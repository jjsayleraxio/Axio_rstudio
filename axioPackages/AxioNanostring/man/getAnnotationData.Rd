\name{getAnnotationData}
\alias{getAnnotationData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Annotation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAnnotationData(x, multiplexed = FALSE)
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
    a <- x$sampleAttributes
    b <- x$laneAttributes
    names(a) <- names(b) <- c("A", "B")
    x <- rbind(a, b)
    cMat <- data.frame(matrix(nrow = length(xNames), ncol = nrow(x)))
    rownames(cMat) <- xNames
    colnames(cMat) <- x[, 1]
    for (i in seq(1, length(xNames))) {
        cMat[i, ] <- x[, 2]
    }
    cMat$Date <- as.Date(cMat$Date, format = "\%Y\%m\%d")
    return(cMat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
