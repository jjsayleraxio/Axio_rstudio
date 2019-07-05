\name{plotGWASPValues}
\alias{plotGWASPValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to create a points plot or hexbin plot of GWAS results across chromosomes}
\description{
Plots -log10 scaled p-values against chromosomal position using points for each p-value or hexagonal bins.  User can also insert lines for significance thresholds and request labels for both the lines and significant SNP IDs.
}
\usage{
plotGWASPValues(position, chromosome, pValue, geneIDs, sort = TRUE, type = c("points", "hexbin"), threshold = NULL, cutoff = NULL, suggestive = NULL, threshLabels = FALSE, geneLabels = FALSE, mainTitle = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{Vector of chromosomal positions.  If the positions are not sorted, then sort must be set to TRUE.}
  \item{chromosome}{Vector of chromosome SNP resides on.}
  \item{pValue}{Vector of -log10 transformed p-values}
  \item{geneIDs}{Vector of SNP IDs.}
  \item{sort}{Logical value undicating whether sorting of chromosomal position should occur.}
  \item{type}{Either points or hexbin.}
  \item{threshold}{Threshold defining significant genes.  Must be specified in raw, untransformed scale.}
  \item{cutoff}{P-value for removing  genes for display.  If observed p-value is below this threshold, the gene is not displayed on the graph.  Must be specified in raw, untransformed scale.}
  \item{suggestive}{Threshold defining suggestive significant genes.  Must be specified in raw, untransformed scale.}
  \item{threshLabels}{Logical.  Print labels on the threshold lines.}
  \item{geneLabels}{Logical.  Print labels on the significant genes.}
  \item{mainTitle}{Character string describing the main title of the graph.}
  \item{\dots}{Other parameters passed to call to plot.}
}
\details{
Creates a plot of -log10 transformed p-value versus chromosomal position.  Thresholds for significance and labeling of significant SNP IDs are also possible.
}
\value{
Creates plot, but returns nothing.
}
\author{David Henderson <DNADave@Merck.Com>}
\seealso{\code{\link{plotGWASPValues}},\code{\link{plotGWASPValues.hexbin}}}
\examples{
data(gwasResults)
plotGWASPValues.default( position = gwasResults$Position , chromosome = gwasResults$chr , pValue = -log10( gwasResults$p.value.additive) , geneIDs = gwasResults$dbSNPid , sort = TRUE , type = "points" , threshold = 2.43e-8, cutoff = NULL , suggestive = 1e-5 )
plotGWASPValues.default( position = gwasResults$Position , chromosome = gwasResults$chr , pValue = -log10( gwasResults$p.value.additive) , geneIDs = gwasResults$dbSNPid , sort = TRUE , type = "hexbin" , threshold = 2.43e-8, cutoff = NULL , suggestive = 1e-5 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GWAS}
\keyword{plot}% __ONLY ONE__ keyword per line
