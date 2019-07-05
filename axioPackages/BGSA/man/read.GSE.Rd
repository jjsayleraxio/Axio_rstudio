\name{read.GSE}
\alias{read.GSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{readGSE}
\description{This function reads the gene expression file along with the sample labels.}
\usage{
read.GSE(expressionFile, labelsFile, outFile = "gse.Rdata", BM.attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"), BM.filters = "affy_hg_u133_plus_2", gmt.file = "c5mf.txt", lower = 5, upper = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionFile}{
expressionFile includes the gene expression values. Each row corresponds to one probe, and each column corresponds to one sample. The first column should have the probe ids. 
}
  \item{labelsFile}{
labelsFile should include one line providing the corresponding group (e.g., case/control) for each sample.
}
  \item{outFile}{
outfile is the name of the output file, which must have .Rdata extension. The output file will include the gene expression data (called data), the labels, the gene names, the set identifers (setInd), and the name of pathways.
}
  \item{BM.attributes}{
BM.attributes specities the attributes for the getBM function of the biomaRt package.
}
  \item{BM.filters}{
BM.filters specities the filters for the getBM function of the biomaRt package.
}
  \item{gmt.file}{
gmt.file, is a text file obtained from the GSEA database; it should include the pathways.
}
  \item{lower}{
The minimum number of genes assigned to a pathway. Pathways that include fewer genes will be excluded from the analysis. 
}
  \item{upper}{
The maximum number of genes assigned to a pathway. Pathways that include more genes will be excluded from the analysis. 
}
}
\details{
Before running this function, make sure there are three text files available in the directory: 1) gene expression data (expressionFile), 2) labels (labelsFile), and 3) list of pathways obtained from the GSEA website (gmt.file). My website includes three sample files that show the required formats. 
}
\value{
The output (outFile) is a .Rdata file that includese gene expression data (called data), the labels, the gene names, the set identifers (setInd), and the name of pathways.
}
\references{
Shahbaba B, Tibshirani R, Shachaf, CM, and Plevritis SK (2011), Bayesian gene set analysis for identifying significant biological pathways, Journal of the Royal Statistical Society, Series C, Volume 60, Issue 4, 541-557.
}
\author{
Babak Shahbaba
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
BGSA
}
\examples{
read.GSE('GSE8401_gene_exp.txt', 'GSE8401_labels.txt')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
