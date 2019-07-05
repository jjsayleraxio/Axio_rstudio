\name{computePowerElstonPGx}
\alias{computePowerElstonPGx}
\title{Function to compute power of rejecting Ho in a pharmacogenetic study.}
\description{
Calculates the power to reject Ho given a Type I error rate and sample size in a
pharmacogenetic study and a list of genetic effects.  Uses the method for normally distributed traits first
described by Elston, et. al. (1999), but as presented in Kelly, Stallard, and
Whittaker (2005).
}
\usage{
computePowerElstonPGx(
       effects ,
       sigma = 1 ,
       p = 0.2 ,
       n = 100 ,
       r = 0.5 ,
       gaction = c( "additive" , "dominant" , "recessive" ) ,
       alpha = 0.05 ,
       beta = 0.8 ,
       type = c( "main" , "interaction" ) )
}
\arguments{
  \item{effects}{A list of genetic effects.  A three element vector for main effects and a n by three matrix for interaction effects.  One entry for each effect examined.}
  \item{sigma}{The common variance associated with either the main genotype effect, or the interaction effect.}
  \item{p}{The minor allele frequency.  May be a vector of frequencies.}
  \item{n}{The proposed sample size.}
  \item{r}{The proportion of the study population in the reference treatment arm.}
  \item{gaction}{The mode of genetic response, either "additive" for a linear response, or "dominance", or "recessive" for a nonlinear response.}
  \item{alpha}{The Type-I error rate.}
  \item{beta}{The Type-II error rate.}
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
A powerObject object.  This contains the genetic effect list, the standard deviation examined, the proposed sample size, the proportion in the treatment arm, the type of test, the power of rejecting Ho given the function parameters and the sample size required to meet the given Type II error.
}
\references{
RC Elston, RM Idury, LR Cardon, and JB Lichter.  1999.  The study of candidate genes in drug trials: sample size considerations.  Statistics in Medicone 18:741-751

PJ Kelly, N Stallard, and JC Whittaker.  2005.  Statistical design and analysis of pharmacogenetic trials.  Statistics in Medicine 24:1495-1508
}
\author{David A. Henderson <\email{David.Henderson@Covance.Com}>}
\note{
Still a work in progress.  API may change without notice.
}
\seealso{\code{\link{getPowerElstonPGxNormal}},\code{\link{getSampleElstonPGxNormal}}}
\examples{
geneticMeans <- list( additive = c( 100 , 110 , 120 ) , dominance = c( 115 , 115 , 120 ) )
geneticVariance <- 30^2
computePowerElstonPGx( geneticMeans , geneticVariance , p = seq( 0.05 , 0.2 , 0.05 ) , n = 100 , r = 0.5 ,
                          gaction = "additive" , alpha = 0.05 , beta = 0.8 , type = "main" )
}
\keyword{power calculation}
\keyword{sample size}
