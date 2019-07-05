createGroups <- function( leaveOut = 1 , n = 10 , nFold = 10 )
{
  groups <- vector( mode = "list" , length = nFold )
  samps <- c( 1:n )
  for ( j in 1:nFold )
  {
    groups[[j]] <- sample( samps , size = ifelse( length( samps ) < leaveOut , length( samps ) , leaveOut ) )
    samps <- samps[which( !( samps %in% groups[[j]] ) )]
  }
  return( groups )
}
