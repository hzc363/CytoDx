#' convert a flowSet to a data frame
#'
#' A function that convert a flowSet to a data frame.
#'
#' @param flowSet A flowSet object
#' @param fcsFiles A vector containing the name of each fcs file included in
#'   flowSet.
#' @param y The clinical outcome each fcs file associated with. Null for testing
#'   data.
#' @return Returns a data frame containing the cytometry data. Cells from
#'   different fcs files are combined into one flow frame. A new column,
#'   xSample, is introduced to indicate the origin of each cell. The data frame
#'   also includes the clinical outcome y.
#' @examples
#' library(flowCore)
#' # Find the table containing fcs file names in CytoDx package
#' path=system.file("extdata",package="CytoDx")
#' # read the table
#' fcs_info = read.csv(file.path(path,"fcs_info.csv"))
#' # Specify the path to the cytometry files
#' fn = file.path(path,fcs_info$fcsName)
#' fSet = read.flowSet(fn)
#' df = set2DF(flowSet=fSet,fcsFiles=fn,y = fcs_info$Label)
#' @export


set2DF=function(flowSet,fcsFiles,y = NULL){
  expr=flowCore::fsApply(flowSet,function(x){
    v=flowCore::exprs(x);
    return(v)
  })
  eventN=flowCore::fsApply(flowSet,function(x){
    v=flowCore::exprs(x);
    n=nrow(v)
    return(n)
  })
  Label=flowCore::fsApply(flowSet,function(x){
    v=flowCore::pData(flowCore::parameters(x));
    return(v)
  })
  #annotate expr colnames
  channels=Label[[1]]$name
  antibodies=Label[[1]]$desc
  #antibodies=toupper(antibodies)
  #channels=toupper(channels)
  antibodies=sapply(1:length(antibodies), function(i){
    if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
  })
  colnames(expr)=antibodies
  #add sample_id
  sample_id=rep(fcsFiles,times = eventN)
  if(is.null(y)){y=rep(NA,nrow(expr))}else{
    y = as.data.frame(y)
    y=y[rep(1:nrow(y),times = eventN),,drop=FALSE]
  }
  expr=cbind.data.frame(expr,"xSample"=sample_id,y)
  return(expr)
}
