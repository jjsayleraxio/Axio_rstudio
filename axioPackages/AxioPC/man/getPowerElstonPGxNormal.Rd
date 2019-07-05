\name{getPowerElstonPGxNormal}
\alias{getPowerElstonPGxNormal}
\title{Function to compute power of rejecting Ho in a pharmacogenetic study.}
\description{
Calculates the power to reject Ho given a Type I error rate and sample size in a
pharmacogenetic study.  Uses the method for normally distributed traits first
described by Elston, et. al. (1999), but as presented in Kelly, Stallard, and
Whittaker (2005).
}
\usage{
getPowerElstonPGxNormal( gmeans,
                         sigma2,
                         p = 0.2,
                         n = 100,
                         r = 0.5,
                         gaction = c("additive", "dominant", "recessive"),
                         alpha = 0.05,
                         type = c("main", "interaction") )
}
\arguments{
  \item{gmeans}{A vector (for main genetic effects) or matrix (for gene by treatment interaction effects) of genotype means.  If a matrix of interaction effects, each row or column corresponds to a single treatment response.}
  \item{sigma2}{The common variance associated with either the main genotype effect, or the interaction effect.}
  \item{p}{The minor allele frequency.}
  \item{n}{The proposed sample size.}
  \item{r}{The proportion of the study population in the reference treatment arm.}
  \item{gaction}{The mode of genetic response, either "additive" for a linear response, or "dominance", or "recessive" for a nonlinear response.}
  \item{alpha}{The Type-I error rate.}
  \item{type}{The type of test performed: either "main" for main genotype effects, or "interaction" for genotype by treatment interaction effects.}
}
\details{
  This function uses the method described in Kelly, Stallard, and Whittaker (2005)
  to calculate the power to reject Ho given a prespecified sample size and Type I
  error rate.  The algorithm first finds Z scores corresponding to the appropriate
  two sided hypothesis test.  The genotypic frequencies are calculated assuming
  Hardy-Weinberg equilibrium.  The distribution of genotypes across treatments
  follows the proportion described by the single parameter \code{r}.  Fixed
  coefficients that describe contrasts corresponding to different types of gene
  action as described in the Kelly, Stallard, and Whittaker paper are also used.
  
  The formula used in the calculations is \deqn{%
    E_A = \frac{\sum{i=0}{2} \sum{j=0}{J} w_{ij}\mu_{ij}}{sqrt{\sigma^2 \sum{i=0}{2} \sum{j=0}{J} frac{w_{ij}^2}{n_{ij}}}}%
    \\ 1 - \beta = \Phi\left( E_A - z_{\alpha/2} \right) + \Phi\left( -E_A - z_{\alpha/2} \right)%
  }{%
    found in the Kelly, Stallard, and Whittaker paper.%
  }
}
\value{
The power of rejecting Ho given the function parameters.
}
\references{
RC Elston, RM Idury, LR Cardon, and JB Lichter.  1999.  The study of candidate genes in drug trials: sample size considerations.  Statistics in Medicone 18:741-751

PJ Kelly, N Stallard, and JC Whittaker.  2005.  Statistical design and analysis of pharmacogenetic trials.  Statistics in Medicine 24:1495-1508
}
\author{David A. Henderson <\email{David.Henderson@Covance.Com}>}
\note{Still a work in progress.  API may change without notice.}
\seealso{\code{\link{getSampleElstonPGxNormal}}}
\examples{
geneticMeans <- matrix( c( 100 , 110 , 120 , 115 , 115 , 120 ) , nrow = 2 , byrow = TRUE )
geneticVariance <- 30^2
getPowerElstonPGxNormal( geneticMeans , geneticVariance , p = 0.3 , n = 100 , r = 0.5 ,
                          gaction = "additive" , alpha = 0.05 , type = "interaction" )
}
\keyword{power calculation}
\keyword{sample size}
