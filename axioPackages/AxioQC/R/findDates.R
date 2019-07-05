findDates <- function( x , type )
{
  if ( type == "quarterly" )
  {
    byVar <- "3 months"
    Year <- list( Q1 = c( "Jan" , "Feb" , "Mar" ) ,
        Q2 = c( "Apr" , "May" , "Jun" ) ,
        Q3 = c( "Jul" , "Aug" , "Sep" ) ,
        Q4 = c( "Oct" , "Nov" , "Dec" ) )
    minD <- format( min( x , na.rm = TRUE ) , "%b" )
    maxD <- format( max( x , na.rm = TRUE ) , "%b" )
    minY <- format( min( x , na.rm = TRUE ) , "%Y" )
    maxY <- format( max( x , na.rm = TRUE ) , "%Y" )
    minQuarter <- which( unlist( lapply( Year , function( x , minD ) minD %in% x , minD ) ) )
    maxQuarter <- which( unlist( lapply( Year , function( x , maxD ) maxD %in% x , maxD ) ) )
    startDate <- as.Date( paste( minY , c( "01" , "04" , "07" , "10" )[minQuarter] , "01" , sep = "-" ) )
    endDate <- as.Date( paste( maxY , c( "03" , "06" , "09" , "12" )[maxQuarter] , c( "31" , "30" , "30" , "31" )[maxQuarter] , sep = "-" ) )
    dateSeq <- seq( startDate , endDate + 1 , byVar )
    qList <- vector( mode = "list" )
    for ( i in seq( 2 , length( dateSeq ) ) )
    {
      bDate <- format( dateSeq[i - 1] , "%b" )
      name <- paste( names( Year[which( unlist( lapply( Year , function( x , bDate ) bDate %in% x , bDate ) ) )] ) , format( dateSeq[i - 1] , "%Y" ) , sep = "-" )
      qList[[name]] <- c( dateSeq[i - 1] , dateSeq[i] - 1 )
    }
    return( qList )
  }
  else if ( type == "bimonthly" )
  {
    byVar <- "2 months"
    Year <- list( B1 = c( "Jan" , "Feb" ) ,
        B2 = c( "Mar" , "Apr" ) ,
        B3 = c( "May" , "Jun" ) ,
        B4 = c( "Jul" , "Aug" ) ,
        B5 = c( "Sep" , "Oct" ) ,
        B6 = c( "Nov" , "Dec" ) )
    minD <- format( min( x , na.rm = TRUE ) , "%b" )
    maxD <- format( max( x , na.rm = TRUE ) , "%b" )
    minY <- format( min( x , na.rm = TRUE ) , "%Y" )
    maxY <- format( max( x , na.rm = TRUE ) , "%Y" )
    minQuarter <- which( unlist( lapply( Year , function( x , minD ) minD %in% x , minD ) ) )
    maxQuarter <- which( unlist( lapply( Year , function( x , maxD ) maxD %in% x , maxD ) ) )
    startDate <- as.Date( paste( minY , c( "01" , "04" , "07" , "10" )[minQuarter] , "01" , sep = "-" ) )
    endDate <- as.Date( paste( maxY , c( "03" , "06" , "09" , "12" )[maxQuarter] , c( "31" , "30" , "30" , "31" )[maxQuarter] , sep = "-" ) )
    dateSeq <- seq( startDate , endDate + 1 , byVar )
    qList <- vector( mode = "list" )
    for ( i in seq( 2 , length( dateSeq ) ) )
    {
      bDate <- format( dateSeq[i - 1] , "%b" )
      name <- paste( names( Year[which( unlist( lapply( Year , function( x , bDate ) bDate %in% x , bDate ) ) )] ) , format( dateSeq[i - 1] , "%Y" ) , sep = "-" )
      qList[[name]] <- c( dateSeq[i - 1] , dateSeq[i] - 1 )
    }
    return( qList )
  }
  else if ( type == "monthly" )
  {
    byVar <- "month"
    startDate <- as.Date( paste( months( min( x ) ) , "-01-" , format( min( x ) , "%Y" ) , sep = "" ) , format = "%b-%d-%Y" )
    endMonth <- switch( months( max( x ) ) , January = "February" , February = "March" , March = "April" , April = "May" , May = "June" , June = "July" , July = "August" , August = "September" , September = "October" , October = "November" , November = "December" , December = "January" )
    endYear <- switch( endMonth , December = paste( as.numeric( format( max( x ) , "%Y" ) ) + 1 ) , format( max( x ) , "%Y" ) )
    endDate <- as.Date( paste( endYear , endMonth , "01" , sep = "-" ) , format = "%Y-%b-%d" ) - 1
    dateSeq <- seq( startDate , endDate , byVar )
    qList <- vector( mode = "list" )
    for ( i in seq( 2 , length( dateSeq ) ) )
    {
      name <- format( dateSeq[i - 1] , "%b-%Y" )
      qList[[name]] <- c( dateSeq[i - 1] , dateSeq[i] - 1 )
    }
    return( qList )
  }
  else if ( type == "biweekly" )
  {
    byVar <- "2 weeks"
    startDate <- switch( weekdays( min( x ) ) , Sunday = min( x - 0 ) , Monday = min( x - 1 ) , Tuesday = min( x - 2 ) , Wednesday = min( x - 3 ) , Thursday = min( x - 4 ) , Friday = min( x - 5 ) , Saturday = min( x - 6 ) )
    endDate <- switch( weekdays( max( x ) ) , Sunday = max( x + 6 ) , Monday = max( x + 5 ) , Tuesday = max( x + 4 ) , Wednesday = max( x + 3 ) , Thursday = max( x + 2 ) , Friday = max( x + 1 ) , Saturday = max( x + 0 ) )
    dateSeq <- seq( startDate , endDate + 1 , byVar )
    qList <- vector( mode = "list" )
    for ( i in seq( 2 , length( dateSeq ) ) )
    {
      bDate <- as.Date( paste( format( dateSeq[i - 1] , "%Y" ) , "01-01" , sep = "-" ) )
      bDate <- switch( weekdays( bDate ) , Sunday = bDate - 0 , Monday = bDate - 1 , Tuesday = bDate - 2 , Wednesday = bDate - 3 , Thursday = bDate + 3 , Friday = bDate + 2 , Saturday = bDate + 1 )
      if ( ( format( dateSeq[i - 1] , "%Y" ) != format( dateSeq[i] , "%Y" ) ) && ( ( dateSeq[i] - as.Date( paste( format( dateSeq[i] , "%Y" ) , "01-01" , sep = "-" ) ) ) > 3 ) )
      {
        name <- paste( "W1" , format( dateSeq[i] , "%Y" ) , sep = "-" )
      }
      else
      {
        week <- paste( "W" , length( seq( bDate , dateSeq[i - 1] , byVar ) ) , sep = "" )
        name <- paste( week , format( dateSeq[i - 1] , "%Y" ) , sep = "-" )
      }
      qList[[name]] <- c( dateSeq[i - 1] , dateSeq[i] - 1 )
    }
    return( qList )
  }
  else if ( type == "weekly" )
  {
    byVar <- "week"
    startDate <- switch( weekdays( min( x ) ) , Sunday = min( x - 0 ) , Monday = min( x - 1 ) , Tuesday = min( x - 2 ) , Wednesday = min( x - 3 ) , Thursday = min( x - 4 ) , Friday = min( x - 5 ) , Saturday = min( x - 6 ) )
    endDate <- switch( weekdays( max( x ) ) , Sunday = max( x + 6 ) , Monday = max( x + 5 ) , Tuesday = max( x + 4 ) , Wednesday = max( x + 3 ) , Thursday = max( x + 2 ) , Friday = max( x + 1 ) , Saturday = max( x + 0 ) )
    dateSeq <- seq( startDate , endDate + 1 , byVar )
    qList <- vector( mode = "list" )
    for ( i in seq( 2 , length( dateSeq ) ) )
    {
      bDate <- as.Date( paste( format( dateSeq[i - 1] , "%Y" ) , "01-01" , sep = "-" ) )
      bDate <- switch( weekdays( bDate ) , Sunday = bDate - 0 , Monday = bDate - 1 , Tuesday = bDate - 2 , Wednesday = bDate - 3 , Thursday = bDate + 3 , Friday = bDate + 2 , Saturday = bDate + 1 )
      if ( ( format( dateSeq[i - 1] , "%Y" ) != format( dateSeq[i] , "%Y" ) ) && ( ( dateSeq[i] - as.Date( paste( format( dateSeq[i] , "%Y" ) , "01-01" , sep = "-" ) ) ) > 3 ) )
      {
        name <- paste( "W1" , format( dateSeq[i] , "%Y" ) , sep = "-" )
      }
      else
      {
        week <- paste( "W" , length( seq( bDate , dateSeq[i - 1] , byVar ) ) , sep = "" )
        name <- paste( week , format( dateSeq[i - 1] , "%Y" ) , sep = "-" )
      }
      qList[[name]] <- c( dateSeq[i - 1] , dateSeq[i] - 1 )
    }
    return( qList )
  }
  else if ( type == "daily" )
  {
    "day"
  }
  else
  {
    stop( "type must be one of: monthly, bimonthly, weekly, quarterly, or daily.\n" )
  }
}

orderDate <- function( x , format = "%b-%d-%Y" )
{
  tmp <- x
  if ( any( unlist( lapply( x , function( x ) grep( "Q" , x , fixed = TRUE ) ) ) ) )
  {
    qs <- c( "Q1" , "Q2" , "Q3" , "Q4" )
    ms <- c( "Jan-01" , "Apr-01" , "Jul-01" , "Oct-01" )
    for ( i in seq( 1 , 4 ) )
    {
      tmp <- gsub( qs[i] , ms[i] , tmp , fixed = TRUE )
    }
  }
  tmp <- as.Date( tmp , format = format )
  return( x[order( tmp )] )
}
