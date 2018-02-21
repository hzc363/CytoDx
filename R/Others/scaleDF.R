#' scale the data
#'
#' A function that performs the rank transformation of the data.
#'
#' @param x A list. Each element should be a numerical matrix, which represents a sample in the data.
#'   Each row of a matrix is a record, each column is a
#'   independent variable included in the model. If your original data is not numerical, using
#'   model.matrix function to convert it to numerical matrices. The names of the list should be the sample names.
#' @return Returns a list containing percentile rank transformed matrix.
#' @examples
#' x = scaleDF(trainData)


scaleDF = function(x,xSample){
  x = cbind.data.frame(x,"xSample"=xSample)%>%
    group_by(xSample)%>%
    mutate_all(.funs=scale)%>%
    as.data.frame()%>%
    select(-xSample)
  return(x)
}
