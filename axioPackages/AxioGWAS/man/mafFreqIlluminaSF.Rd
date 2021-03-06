% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf.R
\name{mafFreqIlluminaSF}
\alias{mafFreqIlluminaSF}
\title{Calculate MAF in parallel}
\usage{
mafFreqIlluminaSF(i, x, Chr = NULL, gender = NULL)
}
\arguments{
\item{i}{index of snps in matrix to calculate MAF for}

\item{x}{matrix of genotypes}

\item{Chr}{vector of which chromosome each single snp is in}

\item{gender}{vector of the genders of the sample for each snp}
}
\value{
vector of named Minor Allele Frequencies
}
\description{
Helper function to fascilitate calculating MAF in parallel
}
\examples{
\dontrun{
 exampleSNP<-matrix(c(0,1,0,1,2,0,0,2,0,0,1,0,1,2,0,0,2,0,0,1,0,1,2,0,0,2,0),nrow=3,byrow=TRUE)

 mafFreqIlluminaSF(c(1:3),exampleSNP,c(1,2,3),c("M","M","F","M","F","F","F","F","M"))
}
}
