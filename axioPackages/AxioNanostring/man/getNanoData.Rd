\name{getNanoData}
\alias{getNanoData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get NanoString Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNanoData(fileCon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileCon}{
%%     ~~Describe \code{fileCon} here~~
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
function (fileCon) 
{
    a <- readLines(fileCon)
    a <- c("<doc>", a, "</doc>")
    b <- xmlParse(a, asText = TRUE)
    r <- xmlRoot(b)
    outPut <- list()
    outPut$header <- read.csv(file = textConnection(xmlValue(r[[1]])), 
        check.names = FALSE)
    outPut$sampleAttributes <- read.csv(file = textConnection(xmlValue(r[[2]])), 
        check.names = FALSE)
    outPut$laneAttributes <- read.csv(file = textConnection(xmlValue(r[[3]])), 
        check.names = FALSE)
    outPut$codeSummary <- read.csv(file = textConnection(xmlValue(r[[4]])), 
        check.names = FALSE)
    return(outPut)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
