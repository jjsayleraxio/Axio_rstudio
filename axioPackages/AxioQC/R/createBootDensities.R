runBootDensities <- function( data , dnarna = "DNA" , date = "Date" , level = "100ng" , type = "quarterly" , nsamples = 20000 , pCutoff = 0.01 , subset = NULL , bootD = TRUE , cumulative = FALSE , LCL = NULL , UCL = NULL )
{
  x <- data[[paste( dnarna , gsub( "ng" , "" , level , fixed = TRUE ) , sep = "" )]]
  if ( !is.null( LCL ) )
  {
    x <- x[-which( x[,level] < LCL ),]
  }
  if ( !is.null( UCL ) )
  {
    x <- x[-which( x[,level] > UCL ),]
  }
  return( createBootDensities( x , date , level , type , nsamples , pCutoff , subset , bootD , cumulative ) )
}

createBootDensities <- function( data , date = "Date" , level = "100ng" , type = "quarterly" , nsamples = 20000 , pCutoff = 0.01 , subset = NULL , bootD = TRUE , cumulative = FALSE )
{
  data <- data[which( !is.na( data[,level] ) ),]
  lower <- pCutoff
  upper <- 1.0 - pCutoff
  boundaries <- findDates( data[,date] , type = type )
  densities <- vector( mode = "list" )
  bootSamples <- vector( mode = "list" )
  if ( bootD )
  {
    bootDensities <- vector( mode = "list" )
  }
  bootPercentiles <- vector( mode = "list" )
  if ( !is.null( subset ) )
  {
    boundaries <- boundaries[subset]
    tmp <- data[0,]
    for ( i in subset )
    {
      tmp <- rbind( tmp , data[data[,date] >= boundaries[[i]][1] & data[,date] <= boundaries[[i]][2],] )
    }
    data <- tmp
  }
  if ( length( boundaries ) > 1 )
  {
    subset <- "All"
    dateNames <- vector()
    if ( cumulative )
    {
      tmp <- data[0,]
    }
    for ( i in names( boundaries ) )
    {
      if ( cumulative )
      {
        tmp <- rbind( tmp , data[which( data[,date] >= boundaries[[i]][1] & data[,date] <= boundaries[[i]][2] ),] )
      }
      else
      {
        tmp <- data[which( data[,date] >= boundaries[[i]][1] & data[,date] <= boundaries[[i]][2] ),]
      }
      if ( nrow( tmp ) > 0 )
      {
        if ( nrow( tmp ) < 2 )
        {
          densities[[i]] <- density( tmp[,level] , na.rm = TRUE , bw = 1 )
        }
        else
        {
          densities[[i]] <- density( tmp[,level] , na.rm = TRUE )
        }
        bootSamples[[i]] <- stack( as.data.frame( smoothedBootstrap( as.vector( na.exclude( tmp[,level] ) ) , nsamples ) ) )$values
        bootPercentiles[[i]] <- findCutoffIntervals( bootSamples[[i]] , c( lower , upper , 0.5 ) )
        if ( bootD )
        {
          bootDensities[[i]] <- density( bootSamples[[i]] , bw = densities[[i]]$bw , n = densities[[i]]$n , na.rm = TRUE )
        }
      }
    }
    densities[["All"]] <- density( data[,level] , na.rm = TRUE )
    bootSamples[["All"]] <- stack( as.data.frame( smoothedBootstrap( as.vector( na.exclude( data[,level] ) ) , nsamples ) ) )$values
    bootPercentiles[["All"]] <- findCutoffIntervals( bootSamples[["All"]] , c( lower , upper , 0.5 ) )
    if ( bootD )
    {
      bootDensities[["All"]] <- density( bootSamples[["All"]] , bw = densities[["All"]]$bw , n = densities[["All"]]$n , na.rm = TRUE )
    }
  }
  else
  {
    if ( is.null( subset ) )
    {
      subset <- "All"
    }
    densities[[subset]] <- density( data[,level] , na.rm = TRUE )
    bootSamples[[subset]] <- stack( as.data.frame( smoothedBootstrap( as.vector( na.exclude( data[,level] ) ) , nsamples ) ) )$values
    bootPercentiles[[subset]] <- findCutoffIntervals( bootSamples[[subset]] , c( lower , upper , 0.5 ) )
    if ( bootD )
    {
      bootDensities[[subset]] <- density( bootSamples[[subset]] , bw = densities[[subset]]$bw , n = densities[[subset]]$n , na.rm = TRUE )
    }
  }
  return( list( boundaries = boundaries , densities = densities , bootSamples = bootSamples , bootPercentiles = bootPercentiles , bootDensities = bootDensities ) )
}

setQuarter <- function( x , type = "quarterly" )
{
  boundaries <- findDates( x , type = type )
  tmp <- vector( length = length( x ) )
  for ( i in names( boundaries ) )
  {
    tmp[which( x >= boundaries[[i]][1] & x <= boundaries[[i]][2] )] <- i
  }
  return( factor( tmp , levels = orderDate( unique( tmp ) ) ) )
}
