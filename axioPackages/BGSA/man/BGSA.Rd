\name{BGSA}
\alias{BGSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{BGSA}
\description{Bayesian Gene Set Analysis}
\usage{
BGSA(file, fout = "outputFile.Rdata", nIter = 1000, burnIn = 200, get.beta = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
This is the name of .Rdata file created by read.GSE
}
  \item{fout}{

}
  \item{nIter}{
Number of MCMC iteration
}
  \item{burnIn}{
Number of burnin samples
}
  \item{get.beta}{
If TRUE, the code provides the posterior samples of regression coefficients.
}
}
\details{
This code uses the .Rdata created by read.GSE() function and performs Bayesian analysis of gene sets.
}
\value{
The output includes the following elements:
post.J: J is a vector of binary latent variables such that J=0 => H_0 is true, J=1 => H_1 is true. Therefore, post.J is the posterior mean of J as a measure for P(H_1 | data)
p.val: This is an alternative measure (comparable to p-value) that we can use to decide whether a gene set is significant
post.beta: this includes the posterior samples for beta's if get.beta was set to TRUE. }
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
read.GSE
}
\examples{
MCMC.res <- BGSA('gse.Rdata', nIter=100, burnIn=20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
