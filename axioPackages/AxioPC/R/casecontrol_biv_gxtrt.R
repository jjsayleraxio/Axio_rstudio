##################
## This code is to generate power for testing the G x TRT trt in Case control study 
## Bivariate outcome
####################
### The model is Logit( P ) = b0 + b1 * TRT + b2 * G + b3 * TRT * G + error
### parameter of interest is b3 = log(ORg x t )
####################
###Assume N00( #non-events in control )
### N01( #events in control ) 
### N10( #non-events in trt )
### N11( #events in trt ) are observable.
### Therefore, the overal risks in trt and control arm are both fixed.
 
### Assume multiplicative (additive on log scale) genetic effect.
####################
### Three prarmeters are provided:

### P: MAF in general population
### ORg x t: Ratio of per-allele OR in the trt arm to that in the control arm. 

### ORg: per-allele OR in control arm. 
# If the gene is assumed that only have interaction effect with the intervention,
# and do NOT have any effect if the patient is in placebo, it could be set as 1.
# sensitivity of the power wrt this parameter can be tested.
###################

## Solve b0, b1, given the risk ratio in control and that in trt arm
###################

getR0 <- function( b0 , b2 , p )
{
    p0 <- ( 1 - p )^2
    p1 <- 2 * p * ( 1 - p )
    p2 <- p * p
    logit <- c( b0 , b0 + b2 , b0 + 2 * b2 ) 
    prob <- exp( logit ) / ( 1 + exp( logit ) )
    return( prob %*% c( p0 , p1 , p2 ) )
}

bisecB0 <- function( a , b , b2 , epslon , p , r0 )
{ 
  while( abs( b - a ) > epslon )
  {
    const <- ( a + b ) / 2
    fa <- getR0( b0 = a , b2 = b2 , p ) - r0
    fc <- getR0( b0 = const , b2 = b2 , p ) - r0
    if ( fa * fc <= 0 )
    {
      b <- const
    }
    if ( fa * fc > 0 )
    {
      a <- const
    }
  }
  return( ( a + b ) / 2 )
}

getR1 <- function( b0 , b2 , b1 , b3 , p )
{
  p0 <- ( 1 - p )^2
  p1 <- 2 * p * (1 - p )
  p2 <- p^2
  logit <- c( b0 + b1 , b0 + b1 + b2 + b3 , b0 + b1 + 2 * ( b2 + b3 ) ) 
  prob <- exp( logit ) / ( 1 + exp( logit ) )
  return( prob %*% c( p0 , p1 , p2 ) )
}

bisecB1 <- function( a , b , b0 , b2 , b3 , epslon , p , r1 )
{ 
  while ( abs( b - a ) > epslon )
  {
    const <- ( a + b ) / 2
    fa <- getR1( b0 , b2 = b2 , b1 = a , b3 = b3 , p ) - r1
    fc <- getR1( b0 , b2 = b2 , b1 = const , b3 = b3 , p ) - r1
    if ( fa * fc <= 0 )
    {
      b <- const
    }
    if ( fa * fc > 0 )
    {
      a <- const
    }
  }
  return( ( a + b ) / 2 )
}

getB0B1 <- function( b2 , b3 , p , r0 , r1 )
{
  b0 <- bisecB0( a = -1000 , b = 1000 , b2 , epslon = 0.0001 , p = p , r0 = r0 )
  b1 <- bisecB1( a = -1000 , b = 1000 , b0 , b2 , b3 , epslon = 0.0001 , p = p , r1 = r1 )
  return( c( b0 , b1 ) )
}

############################
#Inverse Fisher matrix to get se(b3),
#thus get power of detecting non-zero b3
############################

getPowerCC <- function( ORg , ORgxt , N00 , N01 , N10 , N11 , p , alpha )
{
  p0 <- ( 1 - p )^2
  p1 <- 2 * p * ( 1 - p )
  p2 <- p^2
  b2 <- log( ORg )
  b3 <- log( ORgxt )
  N0 <- N00 + N01
  N1 <- N11 + N10
  r0 <- N01 / N0
  r1 <- N11 / N1
  N <- N0 + N1
  b01 <- getB0B1( b2 , b3 , p , r0 , r1 )
  b0 <- b01[1]
  b1 <- b01[2]
  logit <- c( b0 , b0 + b2 , b0 + 2 * b2 , b0 + b1 , b0 + b1 + b2 + b3 , b0 + b1 + 2 * ( b2 + b3 ) ) 
  temp <- exp( -logit ) / ( ( 1 + exp( -logit ) )^2 )
  n.vec <- c( c( p0 , p1 , p2 ) * N0 , c( p0 , p1 , p2 ) * N1 ) 

  fisher <- matrix( NA , ncol = 4 , nrow = 4 )
  fisher[1,1] <- temp %*% n.vec
  fisher[1,2] <- fisher[2,1] <- temp[4:6] %*% n.vec[4:6]
  fisher[1,3] <- fisher[3,1] <- temp[2] * n.vec[2] + temp[5] * n.vec[5] + 2 * 
                 ( temp[3] * n.vec[3] + temp[6] * n.vec[6] )
  fisher[1,4] <- fisher[4,1] <- temp[5] * n.vec[5] + 2 * temp[6] * n.vec[6]
  fisher[2,2] <- temp[4:6] %*% n.vec[4:6]
  fisher[2,3] <- fisher[3,2] <- temp[5] * n.vec[5] + 2 * temp[6] * n.vec[6]
  fisher[2,4] <- fisher[4,2] <- temp[5] * n.vec[5] + 2 * temp[6] * n.vec[6]
  fisher[3,3] <- temp[2] * n.vec[2] + temp[5] * n.vec[5] + 4 * ( temp[3] * n.vec[3] + temp[6] * n.vec[6] )
  fisher[3,4] <- fisher[4,3] <- temp[5] * n.vec[5] + 4 * temp[6] * n.vec[6]
  fisher[4,4] <- temp[5] * n.vec[5] + 4 * temp[6] * n.vec[6]

  se <- sqrt( solve( fisher )[ 4 , 4 ] )
  cutoff <- abs( qnorm( alpha ) ) - b3 / se
  return( 1 - pnorm( cutoff ) )
}

##############################
