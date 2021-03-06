% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{mergeDuplicateData}
\alias{mergeDuplicateData}
\title{Merge Duplicated Data}
\usage{
mergeDuplicateData(x, xNames, uIndx = "PAT_ID", dIndx = "SPCMN")
}
\arguments{
\item{x}{Sample Dataframe, phenotypes}

\item{uIndx}{index of unique IDs}

\item{dIndx}{index of duplicated IDs}

\item{XNames}{the Names of duplicated specimins to merge in dIndx column}

\item{g}{the function to peform the gwas}
}
\value{
data.frame
}
\description{
Merge duplicated values together
}
