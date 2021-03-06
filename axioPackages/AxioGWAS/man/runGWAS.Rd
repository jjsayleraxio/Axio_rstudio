% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas.R
\name{runGWAS}
\alias{runGWAS}
\title{Run a GWAS}
\usage{
runGWAS(pgx, snp, Model, func = gwasFunc, interactionTerm = "ARMACD",
  genderVar = "SEX.Decode", adjust = "BY", replace = NULL,
  snowFall = FALSE, file = NULL, subCoef = NULL, subEffect = NULL,
  dbPath = NULL, dbConnection = NULL, dbTABLENAME = NULL,
  overwriteEntry = NULL, tempFile = tempfile(), ...)
}
\arguments{
\item{pgx}{DF of phenotypes for each subject}

\item{snp}{Matrix of genotype calls, snp by row, subject by column. Genotype format varies by called "func".}

\item{Model}{The original base model for the GWAS, excluding snp and snp interaction terms}

\item{func}{The name of the function to apply in the GWAS. This is meant to be a user-defined function, with the mandatory output of a list containing at least one value called "pval". Inluded in the Package are some examples.}

\item{interactionTerm}{The name of the variable to be included as the interaction term with the snp}

\item{genderVar}{The name of the gender varable in pgx}

\item{adjust}{The method to adjust the P-values by after the GWAS. defaults to "BY", see \code{\link[stats]{p.adjust}} for more options}

\item{replace}{The name of the variable to remove from the Model if necessary}

\item{snowFall}{Boolean Value, Should the GWAS be run in parallel or series. If TRUE, use \code{\link[snowfall]{sfInit}} prior to \code{runGWAS}}

\item{file}{Name of file which to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE}

\item{subEffect}{Effect to keep and to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE}

\item{dbPath}{path to database to save results in if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE}

\item{dbConnection}{a sqlite connection to a DB, alternative to dbPath to save results if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE}

\item{dbTABLENAME}{The name of the table in the DB to save the results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE}

\item{...}{additional arguments passed to func}

\item{subCoeff}{coefficients to keep and to save results to if snp is a  "CNSet" or "SnpSuperSet, and snowFall=TRUE}
}
\value{
a list, with a data.frame of P-Values, and a list of models for each snp
}
\description{
Run the GWAS using the function identified in func. Can be ran in parallel or series
}
\examples{
\dontrun{
examplePGx<-data.frame(var1=c(1,1,1,2,2,2,1,2),var2=c(2,2,2,3,3,3,3,2),var3=c(3,3,3,4,4,4,3,2),gender=c("Male","Female","Male","Female","Male","Female","Male","Male"))
exampleSNP<-matrix(c("A/A","A/B","A/A","A/B","B/B","A/A","A/A","B/B","A/A","A/B","A/A","A/A","B/B","B/B","A/A","B/B"),nrow=2,byrow=TRUE)
exampleSNPInfo<-data.frame(Chr=c(1,2))
exampleModel<-formula("var1~var2+var3+gender")

gwasRes<-runGWAS(pgx=examplePGx,snp=exampleSNP,Model=exampleModel,func = gwasFunc,interactionTerm="var2",genderVar="gender",databasePath="runGWAS.DB", dbTable="runGWAS_Test", snpInfo=exampleSNPInfo, overwriteEntry=TRUE)
}
}
