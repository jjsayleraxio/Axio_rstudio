
makeContrastMatrix <- function( gmeans = NULL , gaction = c( "additive" , "dominance" , "recessive" ) , type = c( "main" , "interaction" ) )
{
  dimIt <- function( x )
  {
    if ( inherits( x , "matrix" ) )
    {
      return( dim( x ) )
    }
    else
    {
      return( c( 1 , length( x ) ) )
    }
  }
  if ( inherits( gmeans , "list" ) )
  {
    dimMat <- matrix( unlist( lapply( gmeans , dimIt ) ) , nrow = length( gmeans ) , byrow = TRUE )
  }
  else
  {
    dimMat <- matrix( dimIt( gmeans ) , nrow = 1 )
  }
  wList <- vector( mode = "list" )
  for ( wi in seq( 1 , nrow( dimMat ) ) )
  {
    if ( type == "main" )
    {
      wList[[wi]] <- switch( gaction ,
          "additive" = matrix( c( -1 , 0 , 1 ) , nrow = 3 ) ,
          "dominance" = matrix( c( -2 , 1 , 1 ) , nrow = 3 ) ,
          "recessive" = matrix( c( -1 , -1 , 2 ) , nrow = 3 ) )
    }
    else if ( type == "interaction" )
    {
      nArms <- dimMat[wi,1] - 1
      if ( gaction == "additive" )
      {
        w <- matrix( NA , nrow = 3 , ncol = nArms + 1 , byrow = TRUE )
        w[,1] <- c( 1 , 0 , -1 )
        for ( i in seq( 2 , nArms + 1 ) )
        {
          w[,i] <- c( -1 , 0 , 1 ) / nArms
        }
        wVec <- as.vector( w[,-1] )
        wVec <- wVec[which( wVec != 0 )]
        wList[[wi]] <- w / min( abs( wVec ) )
      }
      else if ( gaction == "dominance" )
      {
        w <- matrix( NA , nrow = 3 , ncol = nArms + 1 , byrow = TRUE )
        w[,1] <- c( 2 , -1 , -1 )
        for ( i in seq( 2 , nArms + 1 ) )
        {
          w[,i] <- c( -2 , 1 , 1 ) / nArms
        }
        wVec <- as.vector( w[,-1] )
        wVec <- wVec[which( wVec != 0 )]
        wList[[wi]] <- w / min( abs( wVec ) )
      }
      else if ( gaction == "recessive" )
      {
        w <- matrix( NA , nrow = 3 , ncol = nArms + 1 , byrow = TRUE )
        w[,1] <- c( 1 , 1 , -2 )
        for ( i in seq( 2 , nArms + 1 ) )
        {
          w[,i] <- c( -1 , -1 , 2 ) / nArms
        }
        wVec <- as.vector( w[,-1] )
        wVec <- wVec[which( wVec != 0 )]
        wList[[wi]] <- w / min( abs( wVec ) )
      }
      else
      {
        stop( "gaction must be one of: additive, dominance, or recessive" )
      }
    }
    else
    {
      stop( "type must be one of: main or interaction" )
    }
  }
  if ( length( wList ) < 2 )
  {
    wList <- wList[[1]]
  }
  return( wList )
}
