% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas.R
\name{returnDF_writeDB}
\alias{returnDF_writeDB}
\title{Clean GWAS results and store in table}
\usage{
returnDF_writeDB(x, dbPath, TABLE)
}
\arguments{
\item{x}{list of lists generated from gwas_XXX_Funcsf helper functions}

\item{dbPath}{The path to database to save results}

\item{TABLE}{The table in the db to save the results to}
}
\value{
NA
}
\description{
A helper function called by gwas_XXX_FuncSF in runGWAS to extract data and store results in a db
}
