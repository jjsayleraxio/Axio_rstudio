\name{gwasQQPlot}
\alias{gwasQQPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to create a QQ plot of GWAS p-values}
\description{
Draws a QQ plot of random ordered uniform(0,1) values versus observed p-values.
}
\usage{
gwasQQPlot(pValues, main = "QQ Plot", abLine = TRUE, logTransform = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pValues}{Vector or raw, untransformed p-values.}
  \item{main}{Main title of plot.}
  \item{abLine}{Logical.  Draw a line along the 45 degree margin.}
  \item{logTransform}{Logical.  -log10 transform the p-values.}
  \item{...}{Other parameters passed to plot.}
}
\details{
Draws a QQ plot of random ordered uniform(0,1) values versus observed p-values.
}
\value{
Draws plot, returns nothing.
}
\author{David Henderson <DNADave@Merck.Com>}
\seealso{\code{\link{plotGWASPValues}},\code{\link{plotGWASPValues.hexbin}}}
\examples{
data(gwasResults)
gwasQQPlot( pValues = gwasResults$p.value.additive , logTransform = TRUE )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GWAS}
\keyword{plot}% __ONLY ONE__ keyword per line
