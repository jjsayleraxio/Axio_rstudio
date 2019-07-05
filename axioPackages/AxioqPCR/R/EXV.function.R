
EXV.assign <- function( indata , 
			            yvar = "Ct" , 	
			            byvar = "well.set.id" , 
			            y.cutoff = 28.4 , 
			            var.fcn.b = 17.5 ,
			            var.fcn.a.4well = 6.05 * 10^( -14 ) ,
			            var.fcn.a.3well = 2.35 * 10^( -14 ) ,
			            var.fcn.a.2well = 2.23 * 10^( -14 ) ,
			            conc.sd.4well = 0.315 ,
			            conc.sd.3well = 0.122 ,
			            conc.sd.2well = 0.116 ,
			            median.nonNACt = 32 ,
			            dmax.cutoff = 1 )
{
#	names( indata ) <- gsub( as.character( yvar ) , "resp" , names( indata ) )
#	names( indata ) <- gsub( as.character( byvar ) , "byvar" , names( indata ) )

	#----------------------------
	# sample.info;
	#----------------------------
	sample.info <- unique( indata[, -which( names( indata ) == yvar )] )

	# temp.empty is the empty dataset for scenario where there is 0 entry;
	temp.empty <- sample.info
	temp.empty$EXV <- ""
	temp.empty[,paste( yvar , "mean" , sep = "." )] <- NA
	temp.empty$flag.type <- "" 

	#------------------------------
	# check number of NA or 40 Cts;
	#------------------------------
    indata[,yvar] <- as.numeric( indata[,yvar] )
    indata[,yvar][indata[,yvar] == 40] <- NA
    indata$y.is.NA <- as.numeric( is.na( indata[,yvar] ) )
	temp.sum <- DoBy( c( "y.is.NA" , yvar ) , c( byvar ) , data = indata , FUN = c( "sum" , "median" , "mean" , "sd" , "n.get" ) , na.rm = TRUE )

	#=============================================================
	# samples with other than 4 replicates, return error.  No need to follow steps 4 - 10
	#=============================================================
	if ( sum( temp.sum[,paste( yvar , "n.get" , sep = "." )] > 4 | temp.sum$y.is.NA.n.get != 4 ) > 0 )
	{
		print( paste( "Error, some replicates (" , byvar , ") has data with less or more than 4 replicates" , sep = "" ) )
		return( indata[indata[,byvar] %in% c( temp.sum[,byvar][temp.sum[,paste( yvar , "n.get" , sep = "." )] > 4 | temp.sum$y.is.NA.n.get != 4] ),] )
	}

	#==========================================;
	# If 3 Ct are NA, set as EXV and add flag.  Report the remaining Ct.
	#==========================================;
	subset.3.40 <- ( temp.sum$y.is.NA.sum == 3 ) 
	if ( sum( subset.3.40 ) > 0 )
	{
		mean.sd.data.3.40 <- merge( sample.info , temp.sum[subset.3.40,c( byvar , paste( yvar , c( "sd" , "n.get" , "mean" , "median" ) , sep = "." ) )] , by = byvar )
		mean.sd.data.3.40$EXV <- "EXV"
		mean.sd.data.3.40$flag.type <- "3NACt"	
	}
    else
    {
		mean.sd.data.3.40 <- temp.empty[subset.3.40,]
	}

	byvar.exclude.index <- temp.sum[,byvar][subset.3.40]

	#==========================================;
	# If 4 Ct are NA, set as EXV and add flag.  Report 40.
	#==========================================;
	subset.4.40 <- ( temp.sum$y.is.NA.sum == 4 )
	if ( sum( subset.4.40 ) > 0 )
	{
		temp.4.40 <- merge( sample.info , temp.sum[subset.4.40,c( "y.is.NA.sum" , byvar )], by = byvar )
		mean.sd.data.4.40 <- temp.4.40
		mean.sd.data.4.40[,paste( yvar , "mean" , sep = "." )] <- 40
		mean.sd.data.4.40[,paste( yvar , "median" , sep = "." )] <- 40
		mean.sd.data.4.40[,paste( yvar , "sd" , sep = "." )] <- 0
		mean.sd.data.4.40[,paste( yvar , "n.get" , sep = "." )] <- 0
		mean.sd.data.4.40$EXV <- "EXV"
		mean.sd.data.4.40$flag.type <- "4NACt"
		mean.sd.data.4.40 <- mean.sd.data.4.40[,-which( names( mean.sd.data.4.40 ) %in% c( "y.is.NA.sum" , "y.is.NA" ) )]
		byvar.exclude.index <- c( byvar.exclude.index , temp.sum[,byvar][subset.4.40] )
	}
    else
    {
		mean.sd.data.4.40 <- temp.empty[subset.4.40,]
	}

	#==========================================;
	# if 2 Ct are NA, then one can delete NA wells only if the remaining sd
    #  is less than the coryvaronding cutoff.
	#==========================================;
	subset.2.40.can.be.deleted <- ( temp.sum$y.is.NA.sum == 2 ) & 
						          ( ( temp.sum[,paste( yvar , "median" , sep = "." )] <= y.cutoff & temp.sum[,paste( yvar , "sd" , sep = "." )] <= conc.sd.2well ) |
						            ( temp.sum[,paste( yvar , "median" , sep = "." )] > y.cutoff & temp.sum[,paste( yvar , "median" , sep = "." )] <= median.nonNACt & temp.sum[,paste( yvar , "sd" , sep = "." )] <= var.fcn.a.2well * ( temp.sum[,paste( yvar , "median" , sep = "." )] )^( var.fcn.b / 2 ) ) )
	if (sum(subset.2.40.can.be.deleted) > 0) 
	{
		tmp1 <- temp.sum[subset.2.40.can.be.deleted,]
		mean.sd.data.2.40.nonEXV <- merge(sample.info, tmp1[, -which(names(tmp1) %in% c("y.is.NA.sum","y.is.NA.median", "y.is.NA.mean","y.is.NA.sd", "y.is.NA.n.get"))], by = byvar )
		mean.sd.data.2.40.nonEXV$EXV <- ""
		mean.sd.data.2.40.nonEXV$flag.type <- "2NACt, nonEXV"
		byvar.exclude.index <- c( byvar.exclude.index , temp.sum[,byvar][subset.2.40.can.be.deleted] )
	}
    else
    {
		mean.sd.data.2.40.nonEXV <- temp.empty[subset.2.40.can.be.deleted,]
	}

	subset.2.40.EXV <- ( temp.sum$y.is.NA.sum == 2 ) & 
						( ( temp.sum[,paste( yvar , "median" , sep = "." )] <= y.cutoff & temp.sum[,paste( yvar , "sd" , sep = "." )] > conc.sd.2well ) |
						   ( temp.sum[,paste( yvar , "median" , sep = "." )] > y.cutoff & temp.sum[,paste( yvar , "median" , sep = "." )] <= median.nonNACt & temp.sum[,paste( yvar , "sd" , sep = "." )] > var.fcn.a.2well * ( temp.sum[,paste( yvar , "median" , sep = "." )] )^( var.fcn.b / 2 ) ) )
	if ( sum( subset.2.40.EXV ) > 0 ) 
	{
		tmp1 <- temp.sum[subset.2.40.EXV,]
		mean.sd.data.2.40.EXV <- merge( sample.info , tmp1[, -which( names( tmp1 ) %in% c( "y.is.NA.sum" , "y.is.NA.median" , "y.is.NA.mean" , "y.is.NA.sd" , "y.is.NA.n.get" ) )] , by = byvar )
		mean.sd.data.2.40.EXV$EXV <- "EXV"
		mean.sd.data.2.40.EXV$flag.type <- "2NACt"
		byvar.exclude.index <- c( byvar.exclude.index , temp.sum[,byvar][subset.2.40.EXV] )
	}
    else
    {
		mean.sd.data.2.40.EXV <- temp.empty[subset.2.40.EXV,]
	}

	#==========================================;
	# if there is 1 NA Ct;
	# can delete NA wells only if the remaining.sd is less than the coryvaronding cutoff;
	#==========================================;
	subset.1.40.can.be.deleted = (temp.sum$y.is.NA.sum == 1) & 
						((temp.sum$yvar.median <= y.cutoff & temp.sum$yvar.sd <= conc.sd.3well) |
						   (temp.sum$yvar.median > y.cutoff & temp.sum$yvar.median <= median.nonNACt & temp.sum$yvar.sd <= var.fcn.a.3well*(temp.sum$yvar.median)^(var.fcn.b/2)))
	if (sum(subset.1.40.can.be.deleted) > 0) 
	{
		tmp1 = temp.sum[subset.1.40.can.be.deleted, ]
		mean.sd.data.1.40.nonEXV = merge(sample.info, tmp1[, -which(names(tmp1) %in% c("y.is.NA.sum","y.is.NA.median", "y.is.NA.mean","y.is.NA.sd", "y.is.NA.n.get"))], by = byvar )
		mean.sd.data.1.40.nonEXV$EXV = ""
		mean.sd.data.1.40.nonEXV$flag.type= "1NACt, nonEXV"
		byvar.exclude.index = c(byvar.exclude.index, temp.sum[,byvar][subset.1.40.can.be.deleted])
	}
    else
    {
		mean.sd.data.1.40.nonEXV = temp.empty[subset.1.40.can.be.deleted,]
	}

	subset.1.40.EXV = (temp.sum$y.is.NA.sum == 1) & 
						((temp.sum$yvar.median <= y.cutoff & temp.sum$yvar.sd > conc.sd.3well) |
					   (temp.sum$yvar.median > y.cutoff & temp.sum$yvar.median <= median.nonNACt & temp.sum$yvar.sd > var.fcn.a.3well*(temp.sum$yvar.median)^(var.fcn.b/2)))
	if (sum(subset.1.40.EXV) > 0) 
	{
		tmp1 = temp.sum[subset.1.40.EXV,]
		mean.sd.data.1.40.EXV = merge(sample.info, tmp1[, -which(names(tmp1) %in% c("y.is.NA.sum","y.is.NA.median", "y.is.NA.mean","y.is.NA.sd", "y.is.NA.n.get"))], by = byvar )
		mean.sd.data.1.40.EXV$EXV = "EXV"
		mean.sd.data.1.40.EXV$flag.type= "1NACt, EXV"
		byvar.exclude.index = c(byvar.exclude.index, temp.sum[,byvar][subset.1.40.EXV])
	}
    else
    {
		mean.sd.data.1.40.EXV = temp.empty[subset.1.40.EXV,]
	}


	#==========================================;
	# assign NA to 40; calculate median and sd;
	#==========================================;
    indata[,yvar][is.na(indata[,yvar])] = 40
	mean.sd.data = DoBy( yvar , byvar , data = indata[!(indata[,byvar] %in% byvar.exclude.index),] , FUN = c( "median" , "mean" , "max" , "min" , "sd" ) , na.rm = TRUE )

	#=============================
	# first round EXV checking;
	# if CV% > 10%, also EXV;
	#=============================
	mean.sd.data$EXV = ifelse( ( mean.sd.data$yvar.sd / mean.sd.data$yvar.mean ) > 0.1 , "EXV" ,
						ifelse( mean.sd.data$yvar.median <= y.cutoff , 
							ifelse( mean.sd.data$yvar.sd > conc.sd.4well , "EXV" , "" ) ,
								ifelse( mean.sd.data$yvar.sd > var.fcn.a.4well * ( mean.sd.data$yvar.median )^( var.fcn.b / 2 ) , "EXV" , "" ) ) )

	#==========================================;
	# for those are EXV after first round, but with Ct.median < median.nonNACt, check for outlier;
	#=====================================
	# calculate dmax and SD.rem
	#==========================================;
	# for those are not EXV after first round, 
	# or EXV but with Ct.median > median.nonNACt
	# report ;
	# =====================================
	subset.outlier.test = mean.sd.data$EXV == "EXV" & mean.sd.data[,paste( yvar , "median" , sep = "." )] <= median.nonNACt
	if ( sum( !subset.outlier.test ) > 0 )
	{
		mean.sd.data.no.need.check <- merge( sample.info , mean.sd.data[!subset.outlier.test,] , by = byvar )
		mean.sd.data.no.need.check$flag.type <- ifelse( mean.sd.data.no.need.check$EXV == "EXV" , "EXV with median > mean.nonNACt" , "" )
	}
    else
	{
		mean.sd.data.no.need.check <- temp.empty[subset.outlier.test,]
	}

	if ( sum( subset.outlier.test ) > 0 )
	{
		outlier.test.list <- mean.sd.data[subset.outlier.test,]
		outlier.test.data <- merge( outlier.test.list , indata , by = byvar )

		outlier.test.data$max.ct.diff <- ifelse( abs( outlier.test.data[,paste( yvar , "max" , sep = "." )] - outlier.test.data[,paste( yvar , "median" , sep = "." )] ) > abs( outlier.test.data[,paste( yvar , "min" , sep = "." )] - outlier.test.data[,paste( yvar , "median" , sep = "." )] ) , 
								abs( outlier.test.data[,paste( yvar , "max" , sep = "." )] - outlier.test.data[,paste( yvar , "median" , sep = "." )] ), abs( outlier.test.data[,paste( yvar , "min" , sep = "." )] - outlier.test.data[,paste( yvar ,"median" , sep = "." )] ) )
		outlier.test.data$Ct.diff <- outlier.test.data[,yvar] - outlier.test.data[,paste( yvar , "median" , sep = "." )]
		outlier.test.data$to.be.deleted <- ifelse( abs( outlier.test.data$Ct.diff ) == outlier.test.data$max.ct.diff , "ToBeDelete" , "   " )

		mean.sd.remain.input <- outlier.test.data[outlier.test.data$to.be.deleted != "ToBeDelete",]


		if ( length( unique( mean.sd.remain.input[,byvar] ) ) == 1 )
		{
			tmp.y <- mean.sd.remain.input[,yvar][!is.na( mean.sd.remain.input[,yvar] )]
            temp.mean.sd.remain <- data.frame( mean( tmp.y ) , median(tmp.y) ,
                sd( tmp.y ) , length( tmp.y ) , unique( mean.sd.remain.input[,byvar] ) )
            names( temp.mean.sd.remain ) <- c( paste( yvar , "mean" , sep = "." ) , 
                paste( yvar , "median" , sep = "." ) ,
                paste( yvar , "sd" , sep = "." ) ,
                paste( yvar , "n.get" , sep = "." ) ,
                byvar ) 
          }
        else
        {
			temp.mean.sd.remain = DoBy( yvar , byvar , data = mean.sd.remain.input, FUN = c( "mean" , "median" , "sd" , "n.get" ) , na.rm = TRUE )
		}


		mean.sd.remain.data <- merge( temp.mean.sd.remain , outlier.test.data[outlier.test.data$to.be.deleted == "ToBeDelete",c( "max.ct.diff" , byvar )] , by = byvar )
		mean.sd.remain.data$EXV <- ifelse( mean.sd.remain.data[,paste( yvar , "median" , sep = "." )] <= y.cutoff , 
						ifelse( mean.sd.remain.data$max.ct.diff <= dmax.cutoff | mean.sd.remain.data[,paste( yvar , "sd" , sep = "." )] > conc.sd.3well , "EXV" , "" ) ,
						ifelse( mean.sd.remain.data$max.ct.diff <= dmax.cutoff | mean.sd.remain.data[,paste( yvar , "sd" , sep = "." )] > var.fcn.a.3well * ( mean.sd.remain.data[,paste( yvar , "median" , sep = "." )] )^( var.fcn.b / 2 ) , "EXV" , "" ) )

		#---------------------------------------------------------------------------
		# those that fails the outlier test, report the original mean, sd and flag for EXV
		#---------------------------------------------------------------------------
		byvar.index.outlier.not.filtered <- mean.sd.remain.data[,byvar][mean.sd.remain.data$EXV == "EXV"]
		if ( length( byvar.index.outlier.not.filtered ) > 0 )
		{
			mean.sd.data.EXV.not.filtered <- merge( sample.info , mean.sd.data[mean.sd.data[,byvar] %in% byvar.index.outlier.not.filtered,] , by = byvar )
			mean.sd.data.EXV.not.filtered$flag.type <- "EXV, dmax or sd.remain did not pass outlier rule"
		}
        else
		{
			mean.sd.data.EXV.not.filtered <- temp.empty[mean.sd.data[,byvar] %in% byvar.index.outlier.not.filtered,]
		}

		#---------------------------------------------------------------------------
		# those that pass the outlier test, report the mean, sd with distant well deleted
		#---------------------------------------------------------------------------
		byvar.index.outlier.filtered <- mean.sd.remain.data[,byvar][mean.sd.remain.data$EXV == ""]
		if ( length( byvar.index.outlier.filtered ) > 0 )
		{
			mean.sd.data.filtered <- merge( sample.info , mean.sd.remain.data[mean.sd.remain.data[,byvar] %in% byvar.index.outlier.filtered,] , by = byvar )
			mean.sd.data.filtered$flag.type <- "using 3 wells, dmax large enough and sd.remain small enough"
		}
        else
		{
			mean.sd.data.filtered <- temp.empty[mean.sd.remain.data[,byvar] %in% byvar.index.outlier.filtered,]
		}
	}
    else
	{
		mean.sd.data.EXV.not.filtered <- temp.empty[subset.outlier.test,]
		mean.sd.data.filtered <- temp.empty[subset.outlier.test,]
	}

	keep.list <- c( names( sample.info ) , "EXV" , paste( yvar , "mean" , sep = "." ) , "flag.type" )
	mean.sd.data.out <- rbind( mean.sd.data.3.40[,keep.list] ,
                               mean.sd.data.4.40[,keep.list] , 
                               mean.sd.data.2.40.nonEXV[,keep.list] ,
                               mean.sd.data.2.40.EXV[,keep.list] ,
                               mean.sd.data.1.40.nonEXV[,keep.list] ,
                               mean.sd.data.1.40.EXV[,keep.list] ,
                               mean.sd.data.EXV.not.filtered[,keep.list] ,
                               mean.sd.data.filtered[,keep.list] ,
                               mean.sd.data.no.need.check[,keep.list] )


	rm( mean.sd.data.3.40 , mean.sd.data.4.40 , mean.sd.data.2.40.nonEXV ,
        mean.sd.data.2.40.EXV , mean.sd.data.1.40.nonEXV , mean.sd.data.1.40.EXV ,
        mean.sd.data.EXV.not.filtered , mean.sd.data.filtered ,
        mean.sd.data.no.need.check )

#	names(mean.sd.data.out) = gsub("yvar", yvar, names(mean.sd.data.out))
#	names(mean.sd.data.out) = gsub("byvar", byvar, names(mean.sd.data.out))
	return( mean.sd.data.out )
}

#---------------
# customized functions for EXV assignment
#---------------

n.get <- function( x , na.rm = TRUE )
{
  ifelse( ( !na.rm ) , length( x ) , length( x[which( !is.na( x ) )] ) )
} 

#sd <- function( x , na.rm = FALSE ) 
#{ 
#  if ( is.matrix( x ) )
#  {
#    return( apply( x , 2 , my.sd , na.rm = na.rm ) )
#  }
#  else if ( is.vector( x ) )
#  { 
#    if ( na.rm )
#    {
#      x <- x[which( !is.na( x ) )]
#    }
#    if ( length( x ) == 0 )
#    {
#      return( NA )
#    }
#    return( sqrt( var( x , na.rm = na.rm ) ) )
#  } 
#  else if ( is.data.frame( x ) )
#  {
#    return( sapply( x , my.sd , na.rm = na.rm ) )
#  }
#  else
#  { 
#    x <- as.vector( x ) 
#    return( my.sd( x , na.rm = na.rm ) ) 
#  } 
#}
