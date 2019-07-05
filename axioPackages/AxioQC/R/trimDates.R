trimDates <- function( controlData , excludeList )
{
  if ( !is.null( excludeList[["DNA"]] ) )
  {
    if ( !is.null( excludeList[["DNA"]][["DNA10"]] ) )
    {
      for ( i in seq( 1 , nrow( excludeList[["DNA"]][["DNA10"]] ) ) )
      {
        controlData[["DNA10"]] <- controlData[["DNA10"]][which( controlData[["DNA10"]]$QC_Date < excludeList[["DNA"]][["DNA10"]][i,1] | controlData[["DNA10"]]$QC_Date > excludeList[["DNA"]][["DNA10"]][i,2] ),]
      }
    }
    if ( !is.null( excludeList[["DNA"]][["DNA100"]] ) )
    {
      for ( i in seq( 1 , nrow( excludeList[["DNA"]][["DNA100"]] ) ) )
      {
        controlData[["DNA100"]] <- controlData[["DNA100"]][which( controlData[["DNA100"]]$QC_Date < excludeList[["DNA"]][["DNA100"]][i,1] | controlData[["DNA100"]]$QC_Date > excludeList[["DNA"]][["DNA100"]][i,2] ),]
      }
    }
  }
  if ( !is.null( excludeList[["RNA"]] ) )
  {
    if ( !is.null( excludeList[["RNA"]][["RNA10"]] ) )
    {
      for ( i in seq( 1 , nrow( excludeList[["RNA"]][["RNA10"]] ) ) )
      {
        controlData[["RNA10"]] <- controlData[["RNA10"]][which( controlData[["RNA10"]]$QC_Date < excludeList[["RNA"]][["RNA10"]][i,1] | controlData[["RNA10"]]$QC_Date > excludeList[["RNA"]][["RNA10"]][i,2] ),]
      }
    }
    if ( !is.null( excludeList[["RNA"]][["RNA100"]] ) )
    {
      for ( i in seq( 1 , nrow( excludeList[["RNA"]][["RNA100"]] ) ) )
      {
        controlData[["RNA100"]] <- controlData[["RNA100"]][which( controlData[["RNA100"]]$QC_Date < excludeList[["RNA"]][["RNA100"]][i,1] | controlData[["RNA100"]]$QC_Date > excludeList[["RNA"]][["RNA100"]][i,2] ),]
      }
    }
  }
  return( controlData )
}
