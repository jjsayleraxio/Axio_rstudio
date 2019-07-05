# This is the main BGSA function
#' @export
BGSAsnp <- function( beta , setInd , nIter = 2000 , burnIn = 200 )
{
  p <- length( beta )

  nSet <- length( setInd )
  l.s <- NULL
  for ( i in seq( 1 , nSet ) )
  {
    l.s[i] <- length( setInd[[i]] )
  }

  sets <- list()
  G <- NULL
  post.beta <- list()
  post.beta.p <- list()
  count <- 0
  for ( i in seq( 1 , nSet ) )
  {
    indSet <- list()
    ind <- setInd[[i]]
    indSet$ind <- ind[!is.na( ind )]

    indSet$p <- length( indSet$ind )

    indSet$beta <- beta[indSet$ind]

    indSet$tau2 <- 1
    count <- count + 1
    sets[[count]] <- indSet
  }
  nSet <- length(sets)
# J is a vector of binary latent variables such that J=0 => H_0 is true, J=1 => H_1 is true.
#  We can use the posterior mean of J as a measure for P(H_1 | data)
  post.J <- matrix( NA , nrow = nIter , ncol = nSet )
  post.tau2 <- matrix( NA , nrow = nIter , ncol = nSet )
  post.lambda <- rep( NA , nIter )
  post.sc <- matrix( NA , nIter , 2 )
  post.nu <- rep( NA , nIter )
# This is P(H_0|data)
  post.p.H0 <- matrix( NA , nrow = nIter , ncol = nSet )
# This is an alternative measure (comparable to p-value) that we can use to decide whether a gene set is significant
  p.val <- matrix( NA , nrow = nIter , ncol = nSet )
#Initial values for MCMC and fixed parameters
  lambda <- 0.5	# lambda is P(H_1)
  J <- rbinom( nSet , 1 , lambda )
  mu0 <- 0
# lambda ~ Beta(aBeta, bBeta)
  aBeta <- 1
  bBeta <- 1
# tau^2 ~ (1-lambda) x Inv-chi2(nu, sc0)	+ lambda x Inv-chi2(nu, sc0+sc1)
  nu0 <- 5
  sc0 <- 0.1
  sc1 <- 0.2
  muNu <- 1
  sigNu <- 1
  muSc <- 1
  sigSc <- 1
  samp <- list()
#MCMC samples
  for ( iter in seq( 1 , nIter ) )
  {
    print(iter)
    tau2.0 <- NULL
    tau2.1 <- NULL
    for( i in seq( 1 , nSet ) )
    {
      beta.scaled <- sets[[i]]$beta
      nu_n <- nu0 + p
      V <- mean( beta.scaled^2, na.rm = TRUE)
      if ( J[i] == 0 )
      {
        sigma02_n <- ( nu0 * sc0 + p * V ) / ( nu0 + p )
        z <- rchisq( 1 , nu_n )
        tau2 <- nu_n * sigma02_n / z
        tau2.0 <- c( tau2.0 , tau2 )
      }
      else
      {
        sigma02_n <- ( nu0 * ( sc0 + sc1 ) + p * V ) / ( nu0 + p )
        z <- rchisq( 1 , nu_n )
        tau2 <- nu_n * sigma02_n / z
        tau2.1 <- c( tau2.1 , tau2 )
      }
      sets[[i]]$tau2 <- tau2
      post.tau2[iter,i] <- sets[[i]]$tau2
    }
    thisSc0 <- NULL
    if ( is.null( tau2.0 ) )
    {
      thisSc0 <- rgamma( 20 , muSc , sigSc )
      sc0 <- thisSc0[1]
    }
    else
    {
      for ( Rep in seq( 1 , 20 ) )
      {
        sc0 <- getSigTau2( c( tau2.0 , tau2.1 ) , nu0 , sc0 , c( rep( 0 , length( tau2.0 ) ) , rep( sc1 , length( tau2.1 ) ) ) , muSc , sigSc )
        thisSc0[Rep] <- sc0
      }
    }
    thisSc1 <- NULL
    if ( is.null( tau2.1 ) )
    {
      thisSc1 <- rgamma( 20 , muSc , sigSc )
      sc1 <- thisSc1[1]
    }
    else
    {
      for ( rep in seq( 1 , 20 ) )
      {
        sc1 <- getSigTau2( tau2.1 , nu0 , sc1 , sc0 , muSc , sigSc )
        thisSc1[rep] <- sc1
      }
    }
    thisNu0 <- NULL
    for ( rep in seq( 1 , 20 ) )
    {
      nu0 <- getNuTau2( c( tau2.0 , tau2.1 ) , nu0 , c( rep( sc0 , length( tau2.0 ) ) , rep( sc0 + sc1 , length( tau2.1 ) ) ) , muNu , sigNu )
      thisNu0[rep] <- nu0
    }
    post.sc[iter,] <- c( sc0 , sc1 )
    post.nu[iter] <- nu0
    beta.H0 <- rnorm( 1000 , 0 , sqrt( get.richi2( 1000 , nu0 , sc0 ) ) )
    for ( i in seq( 1, nSet ) )
    {
      q0 <- get.dichi2.log( sets[[i]]$tau2 , nu0 , sc0 ) + log( 1 - lambda )
      q1 <- get.dichi2.log( sets[[i]]$tau2 , nu0 , sc0 + sc1 ) + log( lambda )
      m <- apply( cbind( q0 , q1 ) , 1 , max )
      A <- m + log( exp(q0 - m ) + exp( q1 - m ) )
      prob.1 <- exp( q1 - A )
      post.p.H0[iter, i] <- 1 - prob.1
      picked <- rbinom( 1 , 1 , prob.1 )
      if ( picked == 0 )
      {
        J[i] <- 0
      }
      else
      {
        J[i] <- 1
      }
    }
    post.J[iter, ] = J
    if ( iter > 100 )
    {
      newNj1 <- sum( J == 1 )
      lambda <- rbeta( 1 , aBeta + newNj1 , bBeta + nSet - newNj1 )
    }
    post.lambda[iter] <- lambda
    p.val[iter,] <- pchisq( nu0 * sc0 / post.tau2[iter,] , nu0 )
  }
  return( list(
          sets = sets ,
          post.J = colMeans( post.J[burnIn:nIter,] ) ,
          post.p.H0 = colMeans( post.p.H0[burnIn:nIter,] ) ,
          p.val = colMeans( p.val[burnIn:nIter,] )
       )
  )
}
