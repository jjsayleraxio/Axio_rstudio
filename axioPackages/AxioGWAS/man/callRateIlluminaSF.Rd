% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf.R
\name{callRateIlluminaSF}
\alias{callRateIlluminaSF}
\title{Calculate Call rate for Illumina -Helper func}
\usage{
callRateIlluminaSF(i, x, crCutOff)
}
\arguments{
\item{i}{index in x to run call rate test on}

\item{x}{vector of genotypes}

\item{crCutOff}{cutoff for callrate}
}
\value{
vector of snps that failed to pass Call rate test, True means failure
}
\description{
helper function to Calculate call rates of the vector X, return true or false if fails to meets call rate cutoff, via Snowfall
}
\examples{
\dontrun{
  examplePASS<-c(1,2,1,2,NA,3,1,1,3,1,2,1)
  exampleFAIL<-c(NA,NA,NA,NA,NA,NA,1,3,3,1,3,3)

 callRateIllumina(examplePASS,.8)==FALSE
 callRateIllumina(exampleFAIL,.8)==TRUE
}
}
