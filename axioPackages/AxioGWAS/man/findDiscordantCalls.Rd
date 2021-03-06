% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{findDiscordantCalls}
\alias{findDiscordantCalls}
\title{Find Discordant Genotypes}
\usage{
findDiscordantCalls(x, g, uIndx = "PAT_ID", dIndx = "SPCM")
}
\arguments{
\item{x}{Sample Dataframe, phenotypes}

\item{g}{matric of genotypes, snp by row, sample by column}

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
