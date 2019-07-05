getPowerKellyPGxNormal <- function( gmeans , sigma , p = 0.2 , n = 100 , simSamps = 1000 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , type = c( "main" , "interaction" ) )
{
  gCoef <- switch( gaction ,
      additive = c( -1 , 0 , 1 ) ,
      dominance = c( -2 , 1 , 1 ) ,
      recessive = c( -1 , -1 , 2 )
  )
  gfreq <- c( p^2 , 2 * p * ( 1 - p ) , ( 1 - p )^2 )
  if ( type == "interaction" )
  {
    if ( !inherits( gmeans , "matrix" ) )
    {
      gmeans <- as.matrix( gmeans )
    }
    if ( nrow( gmeans ) < 3 )
    {
      gmeans <- t( gmeans )
    }
    if ( length( r ) > 1 )
    {
      if ( ncol( gmeans ) != length( r ) )
      {
        stop( "Dimensions of r and GeneticMeans not conformable" )
      }
    }
    else
    {
      r <- c( r , 1 - r )
    }
    gmeans <- as.vector( gmeans )
    tfreq <- vector( length = length( r ) * length( gfreq ) )
    tcoef <- vector( length = length( r ) * length( gCoef ) )
    for ( i in seq_along( r ) )
    {
      tfreq[( ( i - 1 ) * length( gfreq ) + 1 ):( i * length( gfreq ) )] <- r[i] * gfreq
      tcoef[( ( i - 1 ) * length( gCoef ) + 1 ):( i * length( gCoef ) )] <- gCoef
    }
    gfreq <- tfreq
    gCoef <- tcoef
    nSamp <- round( n * gfreq )
    if ( length( idx <- which( nSamp < 1 ) ) > 0 )
    {
      for ( j in idx )
      {
        nSamp[j] <- 1
        nSamp[j + 2] <- nSamp[j + 2] - 1
      }
    }
    while( sum( nSamp ) < n )
    {
      diffs <- n * gfreq - nSamp
      diffs[which( diffs < 0 )] <- NA
      nSamp[which.max( diffs )] <- nSamp[which.max( diffs )] + 1
    }
    if ( sum( nSamp ) > n )
    {
      nSamp[which.max( nSamp )] <- nSamp[which.max( nSamp )] - 1
    }
    mu <- matrix( rep( gmeans , nSamp ) , ncol = 1 )
    idx <- seq( 1 , length( nSamp ) , 3 )
    trtMat <- matrix( unlist( lapply( c( idx[-1] ) , function( i , nSamp , N ) { a <- rep( 0 , N ); a[nSamp[i]:nSamp[i + 2]] <- 1 ; return( a ) } , cumsum( nSamp ) , n ) ) , nrow = n )
    addVec <- rep( gCoef , nSamp )
    x1 <- matrix( c( rep( 1 , n) , addVec , as.vector( trtMat ) , as.vector( addVec * trtMat ) ) , nrow = n )
    x0 <- matrix( c( rep( 1 , n) , as.vector( trtMat ) ) , nrow = n )
  }
  else
  {
    nSamp <- c( min( ceiling( gfreq[1] * n ) , 1 ) , max( floor( gfreq[2] * n ) , 1 ) , n - sum( c( min( ceiling( gfreq[1] * n ) , 1 ) , max( floor( gfreq[2] * n ) , 1 ) ) ) )
    mu <- matrix( c( rep( gmeans , nSamp ) ) , ncol = 1 )
    x1 <- matrix( c( rep( 1 , n) , rep( gCoef , nSamp ) ) , ncol = 2 )
    x0 <- matrix( c( rep( 1 , n) ) , ncol = 1 )
  }
  nu1 <- ncol( x1 ) - ncol( x0 )
  nu2 <- n - ncol( x1 )
  pow <- 0
  for ( k in seq( 1 , simSamps ) )
  {
    y <- unlist( lapply( seq( 1 , length( nSamp ) ) , function( i , n , g , s ) return( rnorm( n[i] , g[i] , s ) ) , nSamp , gmeans , sigma ) )
    yhat <- x1 %*% solve( t( x1 ) %*% x1 ) %*% t( x1 ) %*% y
    lambda1 <- ( t( mu ) %*% x1 %*% solve( t( x1 ) %*% x1 ) %*% t( x1 ) %*% mu - t( mu ) %*% x0 %*% solve( t( x0 ) %*% x0 ) %*% t( x0 ) %*% mu ) / sigma^2
    lambda2 <- 0
    for ( i in gmeans )
    {
      idx <- which( mu %in% i )
      print( i )
      print( unique( x1[idx,2] ) )
      lambda2 <- lambda2 + unique( x1[idx,2] ) * sum( ( yhat[idx] - i )^2 )
      print( lambda2 )
    }
    lambda2 <- n * lambda2 / sigma^2
    nu2star <- max( c( ( nu2 + lambda2 )^2 / ( nu2 + 2 * lambda2 ) , nu2 ) )
    pow <- pow + pf( qf( 1 - alpha , nu1 , nu2 ) , nu1 , nu2star , lambda1 , lower.tail = FALSE )
  }
  return( round( pow / simSamps , 4 ) )
}

getPowerKellyPGxNormalVec <- function( p = 0.2 , gmeans , sigma2 , n = 100 , simSamps = 1000 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , type = c( "main" , "interaction" ) )
{
  return( as.vector( sapply( gmeans , getPowerKellyPGxNormal , sigma2 , p , n , simSamps , r , gaction , alpha , type ) ) )
}

computePowerKellyPGx <- function ( effects , sigma = 1 , p = 0.2 , n = 100 , simSamps = 1000 , r = 0.5 , gaction = c( "additive" , "dominance" , "recessive" ) , alpha = 0.05 , beta = 0.8 , type = c( "main" , "interaction" ) )
{
  stopifnot( inherits( effects , "list" ) )
  powerMatrix <- matrix( sapply( p , getPowerKellyPGxNormalVec , effects , sigma^2 , n , simSamps , r , gaction , alpha , type ) , nrow = length( p ) , ncol = length( effects ) , dimnames = list( p , names( effects ) ) )
  sampleMatrix <- matrix( sapply( p , getSampleKellyPGxNormalVec , effects , sigma^2 , r , simSamps , gaction , alpha , beta , type ) , nrow = length( p ) , ncol = length( effects ) , dimnames = list( p , names( effects ) ) )
  if ( length( effects ) < 2 )
  {
    powerMatrix <- t( powerMatrix )
    sampleMatrix <- t( sampleMatrix )
  }
  return( new( "powerObject" , effects = effects,
          MAF = p ,
          power = powerMatrix ,
          sampleSize = sampleMatrix ,
          SD = sigma ,
          N = n ,
          alpha = alpha ,
          beta = beta ,
          proportion = r ,
          type = type ,
          stype = "Normal" ) )
}
