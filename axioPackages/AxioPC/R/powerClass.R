setClass( "powerObject" ,
   representation = representation(
   effects = "list" ,
   MAF = "numeric" ,
   power = "matrix" ,
   sampleSize = "matrix" ,
   SD = "numeric" ,
   N = "numeric" ,
   alpha = "numeric" ,
   beta = "numeric" ,
   proportion = "numeric" ,
   type = "character" ,
   stype = "character" ) ,
  prototype = prototype(
    effects = list() , 
    MAF = numeric( 0 ) , 
    power = matrix() , 
    sampleSize = matrix() , 
    SD = numeric( 0 ) , 
    N = numeric( 0 ) , 
    alpha = numeric( 0 ) ,
    beta = numeric( 0 ) ,
    proportion = numeric( 0 ) ,
    type = "" ,
    stype = ""
   ) )

powerObject <- function( effects , p , powerMatrix , sampleMatrix , sigma , n , alpha , beta , r , type , stype )
{
  return( new( "powerObject" ,
   effects = effects ,
   MAF = p ,
   power = t( powerMatrix ) ,
   sampleSize = t( sampleMatrix ) ,
   SD = sigma ,
   N = n ,
   alpha = alpha ,
   beta = beta ,
   proportion = r ,
   type = type ,
   stype = stype ) )
}

print.powerObject <- function( object )
{
  cat( "Study Type:" , object@stype , "\n\n" )
  cat( "Effect sizes examined\n\n" )
  if ( is.null( names( object@effects ) ) )
  {
    effectLevels <- seq( 1 , length( object@effects ) )
  }
  else
  {
    effectLevels <- names( object@effects )
  }
  for ( i in effectLevels )
  {
    cat( i , "\n" )
    tmp <- as.data.frame( object@effects[[i]] )
    if ( ncol( tmp ) < 3 )
    {
      tmp <- t( tmp )
    }
    if ( object@type == "interaction" )
    {
      rownames( tmp ) <- c( "Control" , paste( "Arm" , seq( 1 , nrow( tmp ) - 1 ) ) )
    }
    else
    {
      rownames( tmp ) <- c( "effect" )
    }
    colnames( tmp ) <- c( "AA" , "Aa" , "aa" )
    print( tmp )
    cat( "\n" )
  }
  cat( "Minor allele frequencies examined\n" )
  cat( object@MAF , "\n" )
  if ( object@stype == "Normal" )
  {
    cat( "\nStandard Deviation examined\n" )
    cat( object@SD , "\n" )
  }
  cat( "\nN examined\n" )
  cat( object@N , "\n" )
  cat( "\nType I Error examined\n" )
  cat( object@alpha , "\n" )
  cat( "\nPower examined\n" )
  cat( object@beta , "\n" )
  if ( object@type == "interaction" )
  {
    cat( "\nProportion in treatment arm\n" )
    if ( length( object@proportion ) < 2 )
    {
      tmp <- as.data.frame( matrix( c( object@proportion , 1 - object@proportion ) , nrow = 1 ) )
    }
    else
    {
      tmp <- as.data.frame( matrix( object@proportion , nrow = 1 ) )
    }
    names( tmp ) <- c( "Control" , paste( "Arm" , seq( 1 , ncol( tmp ) - 1 ) ) )
    rownames( tmp ) <- "% in Arm"
    print( tmp )
  }
  cat( "\nType of test\n" )
  cat( object@type , "\n" )
  if ( !is.null( object@power ) )
  {
    cat( "\nPower from calculations\n" )
    if ( ncol( object@power ) < 3 )
    {
      tmp <- t( as.data.frame( object@power ) )
    }
    else
    {
      tmp <- as.data.frame( object@power )
    }
    rownames( tmp ) <- effectLevels
    colnames( tmp ) <- object@MAF
    print( tmp )
  }
  if ( !is.null( object@sampleSize ) )
  {
    cat( "\nSample Size from calculations\n" )
    if ( ncol( object@sampleSize ) < 3 )
    {
      tmp <- t( as.data.frame( object@sampleSize ) )
    }
    else
    {
      tmp <- as.data.frame( object@sampleSize )
    }
    rownames( tmp ) <- effectLevels
    colnames( tmp ) <- object@MAF
    print( tmp )
  }
  invisible( object )
}
setMethod( "show" , "powerObject" , function( object ) print.powerObject( object ) )
