\name{loadNanoString}
\alias{loadNanoString}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load NanoString Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadNanoString(pathRoot = file.path("."), rccFiles = NULL, zipped = FALSE, multiplexed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathRoot}{
%%     ~~Describe \code{pathRoot} here~~
}
  \item{rccFiles}{
%%     ~~Describe \code{rccFiles} here~~
}
  \item{zipped}{
%%     ~~Describe \code{zipped} here~~
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
function (pathRoot = file.path("."), rccFiles = NULL, zipped = FALSE, 
    multiplexed = FALSE) 
{
    if (is.null(rccFiles)) {
        rccFiles <- list.files(pathRoot)
    }
    nanoMat <- matrix(nrow = 0, ncol = 0)
    annoMat <- matrix(nrow = 0, ncol = 0)
    for (rccFile in rccFiles) {
        if (zipped) {
            files <- as.character(unzip(file.path(pathRoot, rccFile), 
                list = TRUE)[, 1])
            fileCon <- unz(file.path(pathRoot, rccFile), files[1])
        }
        else {
            files <- NULL
            fileCon <- file(file.path(pathRoot, rccFile), open = "r")
        }
        datList <- getNanoData(fileCon)
        close(fileCon)
        nanoCounts <- getCountData(datList, multiplexed)
        nanoAnnotation <- getAnnotationData(datList, multiplexed)
        if (!is.null(files)) {
            for (i in seq(2, length(files))) {
                fileCon <- unz(file.path(pathRoot, rccFile), 
                  files[i])
                datList <- getNanoData(fileCon)
                close(fileCon)
                nanoCounts <- cbind(nanoCounts, getCountData(datList, 
                  multiplexed))
                nanoAnnotation <- rbind(nanoAnnotation, getAnnotationData(datList, 
                  multiplexed))
            }
        }
        if (sum(dim(nanoMat)) < 1) {
            nanoMat <- nanoCounts
            annoMat <- nanoAnnotation
        }
        else {
            nanoMat <- cbind(nanoMat, nanoCounts)
            annoMat <- rbind(annoMat, nanoAnnotation)
        }
    }
    protocoldata <- annoMat[, "FovCount", drop = FALSE]
    protocolMetadata <- data.frame(labelDescription = c("Specifiecd FOV count"), 
        row.names = names(protocoldata))
    annoMat <- annoMat[, -which(names(annoMat) \%in\% "FovCount")]
    metaData <- data.frame(labelDescription = c("Owner of Sample", 
        "Comments", "Date", "Filename - .RLF extension", "Filename - .APF extension", 
        "Observed FOV count", "Scanner ID", "Stage position of cartridge lane", 
        "Density of spots in the lane", "Cartridge ID", "Cartridge barcode"), 
        row.names = names(annoMat))
    experimentData <- new("MIAME", name = "CGL", lab = "Covance Genomics Laboratory", 
        contact = "jeannette.nussbaum@covance.com", title = paste(unique(annoMat$GeneRLF)), 
        abstract = "none", other = list(operators = unique(annoMat$Owner), 
            softwareVersion = datList$header[1, 2], fileVersion = names(datList$header)[2]))
    indx <- match(unique(datList$codeSummary[, 2]), datList$codeSummary[, 
        2])
    annotationData <- datList$codeSummary[indx, c(1, 3)]
    annotationData[, 1] <- as.character(annotationData[, 1])
    annotationData[, 2] <- as.character(annotationData[, 2])
    rownames(annotationData) <- datList$codeSummary[indx, 2]
    if (multiplexed) {
        indx <- grep("Endogenous", annotationData$CodeClass, 
            fixed = TRUE)
        annotationData[indx, "CodeClass"] <- substr(annotationData[indx, 
            "CodeClass"], 1, nchar(annotationData[indx, "CodeClass"]) - 
            1)
        indx <- grep("Housekeeping", annotationData$CodeClass, 
            fixed = TRUE)
        annotationData[indx, "CodeClass"] <- substr(annotationData[indx, 
            "CodeClass"], 1, nchar(annotationData[indx, "CodeClass"]) - 
            1)
    }
    annotationData <- new("AnnotatedDataFrame", data = annotationData, 
        varMetadata = data.frame(labelDescription = c("Class of feature", 
            "Accession Number"), row.names = c("CodeClass", "Accession")))
    nanoSet <- ExpressionSet(assayData = nanoMat, phenoData = new("AnnotatedDataFrame", 
        data = annoMat, varMetadata = metaData), protocolData = new("AnnotatedDataFrame", 
        data = protocoldata, varMetadata = protocolMetadata), 
        experimentData = experimentData, featureData = annotationData)
    return(nanoSet)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
