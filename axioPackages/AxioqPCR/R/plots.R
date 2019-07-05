#

is.shingle <- function( x )
{
  return( inherits( x , "shingle" ) )
}

as.factorOrShingle <- function( x , subset = TRUE , drop = FALSE ) 
{
  x <- if (is.numeric(x))
  {
    as.shingle(x)
  }
  else
  {
    as.factor(x)
  }
  return( x[subset, drop = drop] )
}

parseFormula <- function (model, data, dimension = 2, subset = TRUE, groups = NULL, 
    multiple = FALSE, outer = FALSE, subscripts = FALSE, drop = NULL) 
{
  expr2char <- function(x) paste(deparse(x), collapse = "")
  if (inherits(groups, "formula")) {
    groupVar <- as.character(groups)[2]
    groups <- eval(parse(text = groupVar), data, environment(groups))
  }
  if (is.null(drop)) 
    drop <- TRUE
  if (is.list(drop)) {
    drop.unused.cond <- if (is.null(drop$cond)) 
          TRUE
        else drop$cond
    drop.unused.data <- if (is.null(drop$data)) 
          TRUE
        else drop$data
  }
  else {
    drop.unused.cond <- drop
    drop.unused.data <- drop
  }
  parseSide <- function(model) {
    model.vars <- list()
    while (length(model) == 3 && model[[1]] == as.name("+")) {
      model.vars <- c(model.vars, model[[3]])
      model <- model[[2]]
    }
    rev(c(model.vars, model))
  }
  parseCond <- function(model) {
    model <- substitute(~m, list(m = model))[[2]]
    model.vars <- list()
    while (length(model) == 3 && (model[[1]] == as.name("*") || 
          model[[1]] == as.name("+"))) {
      model.vars <- c(model.vars, model[[3]])
      model <- model[[2]]
    }
    rev(c(model.vars, model))
  }
  lrep <- function(x, n) {
    save.attr <- attributes(x)
    x <- rep(x, n)
    attributes(x) <- save.attr
    x
  }
  concat <- function(arglist) {
    if (length(arglist) == 1) 
      arglist[[1]]
    else if (any(sapply(arglist, is.factor))) {
      factor(unlist(lapply(arglist, as.character)))
    }
    else if (any(sapply(arglist, is.shingle))) {
      stop("shingles can not be concatenated")
    }
    else do.call("c", arglist)
  }
  if (!inherits(model, "formula")) 
    stop("model must be a formula object")
  if (multiple && !outer && !is.null(groups)) {
    multiple <- FALSE
    warning("'multiple=TRUE' ignored ('groups' non-null with 'outer=FALSE')")
  }
  ans <- if (dimension == 2) {
        list(left = NULL, right = NULL, condition = NULL, left.name = character(0), 
            right.name = character(0))
      }
      else if (dimension == 3) {
        list(left = NULL, right.x = NULL, right.y = NULL, condition = NULL, 
            left.name = character(0), right.x.name = character(0), 
            right.y.name = character(0))
      }
      else stop(gettextf("invalid dimension '%s'", as.character(dimension)))
  if (length(model) == 3) {
    if (multiple) {
      varsLHS <- parseSide(model[[2]])
      nLHS <- length(varsLHS)
    }
    else {
      varsLHS <- list(model[[2]])
      nLHS <- 1
    }
  }
  else {
    nLHS <- 1
  }
  modelRHS <- model[[length(model)]]
  if (length(modelRHS) == 3 && modelRHS[[1]] == as.name("|")) 
    modelRHS <- modelRHS[[2]]
  env <- environment(model)
  modelRHS <- model[[length(model)]]
  if (length(modelRHS) == 3 && modelRHS[[1]] == as.name("|")) {
    modelRHS.vars <- parseCond(modelRHS[[3]])
    modelRHS <- modelRHS[[2]]
    if (multiple && dimension == 2) {
      varsRHS <- parseSide(modelRHS)
      nRHS <- length(varsRHS)
    }
    else {
      varsRHS <- list(modelRHS)
      nRHS <- 1
    }
    ans$condition <- vector("list", length(modelRHS.vars))
    names(ans$condition) <- sapply(modelRHS.vars, expr2char)
    for (i in seq_along(modelRHS.vars)) {
      ans$condition[[i]] <- lrep(as.factorOrShingle(eval(modelRHS.vars[[i]], 
                  data, env), subset, drop = drop.unused.cond), 
          nLHS * nRHS)
    }
  }
  else if (multiple && dimension == 2) {
    varsRHS <- parseSide(modelRHS)
    nRHS <- length(varsRHS)
  }
  else {
    varsRHS <- list(modelRHS)
    nRHS <- 1
  }
  if (length(model) == 3) {
    ans$left.name <- expr2char(model[[2]])
    ans$left <- lrep(concat(lapply(varsLHS, function(i) {
                  tmp <- eval(i, data, env)
                  if (!is.matrix(tmp)) 
                    tmp <- if (is.factor(tmp) || is.shingle(tmp)) 
                          tmp[subset, drop = drop.unused.data]
                        else tmp[subset]
                  if (inherits(tmp, "POSIXt")) 
                    tmp <- as.POSIXct(tmp)
                  tmp
                })), nRHS)
  }
  if (dimension == 2) {
    tmp <- eval(varsRHS[[1]], data, env)
    if (is.matrix(tmp)) 
      tmp <- as.data.frame(tmp)
    nobs <- if (is.data.frame(tmp)) 
          nrow(tmp)
        else length(tmp)
    if (nLHS == 1 && nRHS == 1) {
      if (is.data.frame(tmp)) 
        ans$right <- tmp[subset, ]
      else ans$right <- if (is.factor(tmp) || is.shingle(tmp)) 
              tmp[subset, drop = drop.unused.data]
            else tmp[subset]
    }
    else {
      ans$right <- concat(lapply(varsRHS, function(i) {
                tmp <- eval(i, data, env)
                tmp <- if (is.factor(tmp) || is.shingle(tmp)) 
                      tmp[subset, drop = drop.unused.data]
                    else tmp[subset]
                tmp <- lrep(tmp, nLHS)
                if (inherits(tmp, "POSIXt")) 
                  tmp <- as.POSIXct(tmp)
                tmp
              }))
    }
    ans$right.name <- expr2char(modelRHS)
    nRows <- length(ans$right)/(nLHS * nRHS)
  }
  else if (dimension == 3 && length(modelRHS) == 3 && (modelRHS[[1]] == 
        "*" || modelRHS[[1]] == "+")) {
    tmp <- eval(modelRHS[[2]], data, env)
    nobs <- length(tmp)
    if (!is.matrix(tmp)) 
      tmp <- if (is.factor(tmp) || is.shingle(tmp)) 
            tmp[subset, drop = drop.unused.data]
          else tmp[subset]
    ans$right.x <- lrep(tmp, nLHS)
    if (inherits(ans$right.x, "POSIXt")) 
      ans$right.x <- as.POSIXct(ans$right.x)
    tmp <- eval(modelRHS[[3]], data, env)
    if (!is.matrix(tmp)) 
      tmp <- if (is.factor(tmp) || is.shingle(tmp)) 
            tmp[subset, drop = drop.unused.data]
          else tmp[subset]
    ans$right.y <- lrep(tmp, nLHS)
    if (inherits(ans$right.y, "POSIXt")) 
      ans$right.y <- as.POSIXct(ans$right.y)
    ans$right.x.name <- expr2char(modelRHS[[2]])
    ans$right.y.name <- expr2char(modelRHS[[3]])
    nRows <- length(ans$right.x)/nLHS
  }
  else stop("invalid model")
  if (nLHS > 1) 
    LHSgroups <- rep(gl(nLHS, nRows, labels = sapply(varsLHS, 
                expr2char)), nRHS)
  if (nRHS > 1) 
    RHSgroups <- gl(nRHS, nRows * nLHS, labels = sapply(varsRHS, 
            expr2char))
  newFactor <- if (nLHS > 1 && nRHS > 1) {
        interaction2(LHSgroups, RHSgroups, sep = lattice.getOption("interaction.sep"))
      }
      else if (nLHS > 1) 
        LHSgroups
      else if (nRHS > 1) 
        RHSgroups
      else NULL
  if (nLHS == 1 && nRHS == 1) {
    if (!is.null(groups)) 
      ans$groups <- groups
    if (subscripts) 
      ans$subscr <- seq_len(nobs)[subset]
  }
  else if (outer) {
    if (!is.null(groups)) 
      ans$groups <- rep(groups, nLHS * nRHS)
    if (!is.null(newFactor)) {
      if (is.null(ans$cond)) 
        ans$condition <- list(newFactor)
      else ans$condition[[length(ans$condition) + 1]] <- newFactor
    }
    else stop("newFactor cannot be NULL; you have found a bug!")
    if (subscripts) 
      ans$subscr <- as.vector(matrix(seq_len(nobs * nLHS * 
                      nRHS), nrow = nobs)[subset, ])
  }
  else {
    if (is.null(groups) && !is.null(newFactor)) 
      ans$groups <- newFactor
    else stop("newFactor != NULL && groups == NULL does not hold; you have found a bug!")
    if (subscripts) 
      ans$subscr <- seq_len(length(newFactor))
    if (length(newFactor) != nRows * nLHS * nRHS) 
      stop("Length check mismatch; you have found a bug!")
  }
  return( ans )
}

setupQPCRVPTree <- function( layout = NULL , nGroups = 1 )
{
  nCellsPerPage <- unlist( layout[1] ) * unlist( layout[2] )
  nPages <- 1
  remainder <- nGroups
  vtlist <- vector( mode = "list" , length = 1 )
  while ( remainder > 0 )
  {
    vplist <- vector( mode = "list" , length = nCellsPerPage )
    for ( i in seq( unlist( layout[1] ) ) )
    {
      for ( j in seq( unlist( layout[2] ) ) )
      {
        vplist[[(i - 1) * unlist( layout[2] ) + j]] <- viewport(
            layout.pos.row = i , layout.pos.col = j , name = cellName( i , j ) )
      }
    }
    vtlist[[nPages]] <- vpTree( viewport( layout = layout ,
            name = cellGridName( nPages ) ) ,
            do.call( "vpList" , vplist ) )
    nPages <- nPages + 1
    remainder <- remainder - nCellsPerPage
  }
  return( vtlist )
}

plotQPCR.data.vp <- function( x , y )
{
  return(
      viewport(
          name = "dataregion" ,
          x = unit( 2 , "lines" ) ,
          y = unit( 1.75 , "lines" ) ,
          width = unit( 1 , "npc" ) - unit( 2.5 , "lines" ) ,
          height = unit( 1 , "npc" ) - unit( 2.75 , "lines" ) ,
          just = c( "left" , "bottom" ) ,
          xscale = range( x , na.rm = TRUE ) + c( -0.05 , 0.05 ) * diff( range( x , na.rm = TRUE ) ) ,
          yscale = range( y , na.rm = TRUE ) + c( -0.05 , 0.05 ) * diff( range( y , na.rm = TRUE ) )
          )
      )
}

removeNAGroups <- function( levels , groups , y )
{
  levels <- unique( levels )
  groups <- as.character( unlist( groups ) )
  for ( lev in levels )
  {
    idx <- which( groups == lev )
    if ( all( is.na( y[idx] ) ) )
    {
      levels <- levels[-which( levels == lev )]
    }
  }
  return( levels )
}

qpcrPlotTitle <- function( title , titleFontsize = 16 )
{
  return( 
      textGrob(
          title ,
          name = "title" ,
          y = unit( 1 , "npc" ) + unit( 0.75 , "lines" ) ,
          gp = gpar( fontsize = titleFontsize ) ,
          vp = "dataregion"
      )
  )
}

cellName <- function( i , j )
{
  return( paste( "cell" , paste( i , j , sep = "." ) , sep = "-" ) )
}

cellGridName <- function( nPages )
{
  return( paste( "cellgrid" , nPages , sep = "-" ) )
}

cellPath <- function( i , j , cellgrid = NULL )
{
  return( vpPath( cellgrid , paste( "cell" , paste( i , j , sep = "." ) , sep = "-" ) ) )
}

plotQPCRPane <- function( x , y , title , name = NULL , draw = TRUE ,
                          gp = gpar() , yrange = c( 0 , 1 ) , vp = NULL ,
                          xFontsize = 14 , yFontsize = 14 , titleFontsize = 16 ,
                          LLOQ = 34.10 , xAxisVals = NULL , xlab = "X Axis" ,
                          ylab = "Y Axis" )
{
  Colors <- rep( "blue" , length( y ) )
  Colors[which( y > LLOQ )] <- "red"
  qpcrPlot <- gTree(
      x = x ,
      y = y ,
      title = title ,
      name = name ,
      childrenvp = plotQPCR.data.vp( xAxisVals[,2] , yrange ) ,
      children = gList(
          rectGrob( name = "border" , vp = "dataregion" ) ,
          xaxisGrob( name = "xaxis" , vp = "dataregion" , at = xAxisVals[,2] ,
              label = xAxisVals[,1] , gp = gpar( fontsize = round( 0.6 * xFontsize ) ) ) ,
          yaxisGrob( name = "yaxis" , vp = "dataregion" ,
              gp = gpar( fontsize = round( 0.6 * yFontsize ) ) ) ,
          pointsGrob( x , y , name = "points" , vp = "dataregion" ,
              pch = 2 , size = unit( titleFontsize / 28 , "char" ) ,
              gp = gpar( col = Colors ) ) ,
          textGrob( xlab , y = unit( -2 , "lines" ) , name = "xlab" ,
              gp = gpar( fontsize = xFontsize ) , vp = "dataregion" ) ,
          linesGrob( x = unit( range( xAxisVals[,2] , na.rm = TRUE ) + c( -0.05 , 0.05 )
              * diff( range( xAxisVals[,2] , na.rm = TRUE ) ) , "native" ) ,
              y = unit( rep( LLOQ , 2 ) , "native" ) , gp = gpar( col = "red" ) ,
              default.units = "native" , name = "cutoff" , vp = "dataregion" ) ,
          textGrob( ylab , x = unit( -2 , "lines" ) , name = "ylab" ,
              gp = gpar( fontsize = xFontsize ) , rot = 90 , vp = "dataregion" ) ,
          qpcrPlotTitle( title , titleFontsize )
      ) ,
      gp = gp ,
      vp = vp ,
      cl = "plotQPCR"
  )
  if ( draw )
  {
    grid.draw( qpcrPlot )
  }
  return( qpcrPlot )
}

setupQPCRVPPage <- function( nPages = 1 )
{
  nCellsPerPage <- 2
  npages <- 1
  vtlist <- vector( mode = "list" , length = nPages )
  while ( nPages > 0 )
  {
    vplist <- vector( mode = "list" , length = nCellsPerPage )
    vplist[[1]] <- viewport( layout.pos.row = 1 , layout.pos.col = 1 ,
        y = unit( 1 , "npc" ) - unit( 0.75 , "lines" ) ,
        height = unit( 2 , "lines" ) ,
        name = paste( "pageTitle" , npages , sep = "-" ) )
    vplist[[2]] <- viewport( layout.pos.row = 2 , layout.pos.col = 1 ,
        height = unit( 0.95 , "npc" ) ,
        y = unit( 0.95 , "npc" ) - unit( 0.025 , "lines" ) ,
        just = "top" , name = paste( "plotGrid" , npages , sep = "-" ) )
    vtlist[[npages]] <- vpTree( viewport( layout = grid.layout( 2 , 1 ) ,
            name = paste( "pagePlot" , npages , sep = "-" ) ) ,
        do.call( "vpList" , vplist ) )
    npages <- npages + 1
    nPages <- nPages - 1
  }
  return( vtlist )
}

plotQPCR <- function( Data = NULL , Formula = Ct ~ Sample.Mass | Detector.Name ,
    group = "Date" , ylab = "Ct", lty = 3 , xlab = "Sample Mass", LLOQ = 34.10 ,
    main = "UHR Ct vs. Mass Plot: by Date" , col = "Grey" , draw = TRUE ,
    xAxisVal = c( 0.0005 , 0.005 , 0.05 , 0.5 , 5 , 50 ) , layout = c( 3 , 3 ) ,
    xFontsize = 14 , yFontsize = 14 , titleFontsize = 16 , setVP = FALSE )
{
  plotData <- parseFormula( Formula , data = Data , groups = group )
  if ( !is.null( plotData$condition ) )
  {
    grid.newpage()
    layout <- grid.layout( layout[1] , layout[2] )
    groups <- removeNAGroups( levels( unlist( plotData$condition ) ) , plotData$condition , plotData$left )
    nGroups <- length( groups )
    nCellsPerPage <- unlist( layout[1] ) * unlist( layout[2] )
    nPages <- 1
    remainder <- nGroups
    plotVPTree <- setupQPCRVPTree( layout = layout , nGroups = nGroups )
    yrange <- range( plotData$left , na.rm = TRUE )
    xLevels <- unique( plotData$right )
    xLevels <- xLevels[order( as.numeric( xLevels ) )]
    txAxisVal <- xAxisVal[xAxisVal %in% xLevels]
    txAxisVal <- cbind( txAxisVal , seq( 1 , length( txAxisVal ) ) )
    gtree <- vector( mode = "list" , length = 1 )
    while ( remainder > 1 )
    {
      glist <- vector( mode = "list" , length = min( remainder , nCellsPerPage ) )
      pushViewport( plotVPTree[[nPages]][[2]] )
      for ( i in seq( min( ceiling( remainder / unlist( layout[2] ) ) ,
              unlist( layout[1] ) ) ) )
      {
        for ( j in seq( min( remainder , unlist( layout[2] ) ) ) )
        {
          grp <- groups[nGroups - ( remainder - 1 )]
          idx <- unlist( plotData$condition ) == grp
          x <- plotData$right[idx]
          for ( xl in 1:nrow( txAxisVal ) )
          {
            x[x == txAxisVal[xl,1]] <- txAxisVal[xl,2]
          }
          y <- plotData$left[idx]
          glist[[(i - 1) * unlist( layout[2] ) + j]] <- plotQPCRPane( x , y ,
              title = grp , name = paste( "plot" , paste( i , j , sep = "." ) ,
              sep = "-" ) , draw = FALSE , gp = gpar() , vp = cellPath( i , j ,
              cellGridName( nPages ) ) , yrange = yrange , xFontsize = xFontsize ,
              yFontsize = yFontsize , titleFontsize = titleFontsize ,
              xAxisVals = txAxisVal , xlab = xlab , ylab = ylab , LLOQ = LLOQ )
          remainder <- remainder - 1
        }
      }
      if ( setVP )
      {
        plotVP <- viewport( layout.pos.row = 2 , layout.pos.col = 1 ,
            height = unit( 0.95 , "npc" ) ,
            y = unit( 0.95 , "npc" ) - unit( 0.025 , "lines" ) ,
            just = "top" , name = paste( "plotGrid" , nPages , sep = "-" ) )
        gtree[[nPages]] <- gTree( name = paste( "plotGrid" , nPages ,
                sep = "-" ) , childrenvp = plotVPTree[[nPages]] ,
            children = do.call( "gList" , glist ) ,
            layout.pos.col = 1 , layout.pos.row = 1 ,
            vp = plotVP )
      }
      else
      {
        gtree[[nPages]] <- gTree( name = paste( "plotgrid" , nPages ,
                sep = "-" ) , childrenvp = plotVPTree[[nPages]] ,
            children = do.call( "gList" , glist ) )
      }
      nPages <- nPages + 1
    }
  }
  if ( draw )
  {
    for ( i in 1:length( gtree ) )
    {
      if ( draw )
      {
        grid.newpage()
      }
      grid.draw( gtree[[i]] )
    }
  }
  invisible( gtree )
}

plotCt <- function( Data = NULL , Formula = Ct ~ Sample.Mass | Detector.Name ,
    group = "Date" , ylab = "Ct", lty = 3 , xlab = "Sample Mass", LLOQ = 34.10 ,
    main = "UHR Ct vs. Mass Plot: by Date" , col = "Grey" , draw = TRUE ,
    xAxisVal = c( 0.0005 , 0.005 , 0.05 , 0.5 , 5 , 50 ) , layout = c( 3 , 3 ) ,
    xFontsize = 8 , yFontsize = 8 , titleFontsize = 10 )
{
  plots <- plotQPCR( Data = Data , xFontsize = xFontsize , yFontsize = yFontsize ,
      titleFontsize = titleFontsize , Formula = Formula , group = group ,
      xlab = xlab , ylab = ylab , lty = lty , LLOQ = LLOQ , main = main ,
      col = col , draw = FALSE , xAxisVal = xAxisVal , layout = layout ,
      setVP = TRUE )
  pageVPTree <- setupQPCRVPPage( nPages = length( plots ) )
  plotlist <- vector( mode = "list" , length = length( plots ) )
  for ( i in 1:length( plots ) )
  {
    vplist <- vector( mode = "list" , length = 2 )
    textVP <- viewport( name = paste( "pageTitle" , i , sep = "-" ) ,
        y = unit( 1 , "npc" ) - unit( 0.75 , "lines" ) ,
        layout.pos.col = 1 , layout.pos.row = 1 )
    vplist[[1]] <- textGrob( main , name = paste( "pageTitle" , i , sep = "-" ) ,
        vp = textVP )
    vplist[[2]] <- textGrob( unique( Data[,group] ) , name = paste( "pageDate" , i , sep = "-" ) ,
        gp = gpar( fontsize = 0.5 * titleFontsize ) , x = unit( 0.9 , "npc" ) ,
        just = "right" , vp = textVP )
    vplist[[3]] <- plots[[i]]
    plotlist[[i]] <- gTree( name = paste( "pagePlot" , i ,
            sep = "-" ) , childrenvp = pageVPTree[[i]] ,
        children = do.call( "gList" , vplist ) )
  }
  if ( draw )
  {
    for ( i in 1:length( plotlist ) )
    {
      grid.newpage()
      grid.draw( plotlist[[i]] )
    }
  }
  invisible( plotlist )
}

#a <- plotCt( Data = subset( data , Sample.Type == "UHR" ) )
