% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf.R
\name{ifelse1}
\alias{ifelse1}
\title{Alternative ifelse}
\usage{
ifelse1(x, a, b)
}
\arguments{
\item{x}{boolean statement}

\item{a}{outcome if x is TRUE}

\item{b}{outcome if x is FALSE}
}
\value{
value of a or b
}
\description{
Custom ifelse statement
}
\examples{
\dontrun{
ifelse1(TRUE,1,0)==1
ifelse1(FALSE,1,0)==0
}
}
