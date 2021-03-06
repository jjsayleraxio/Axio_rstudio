% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plink.R
\name{plinkifyGenotypes}
\alias{plinkifyGenotypes}
\title{Plink-ify genes}
\usage{
plinkifyGenotypes(x)
}
\arguments{
\item{x}{numeric genotype}
}
\description{
convert numeric genotypes to plink format
}
\examples{
\dontrun{
 exampleSNP<-c(1,2,1,2,3,1,1,3,1)
 
 plinkifyGenotypes(exampleSNP[1])=="A A"  
}
}
