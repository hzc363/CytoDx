#' Percentile rank transformation of the data
#'
#' A function that performs the rank transformation of the data.
#'
#' @param x A data frame containing the pooled data from fcs files. Each row is
#'   a cell, each column is a marker.
#' @param xSample A vector specifying which sample each cell belongs to. Length
#'   must equal to nrow(x).
#' @return Returns data frame containing rank transformed data.
#' @importFrom dplyr mutate_all
#' @examples
#' x <- pRank(x=iris[,1:4],xSample=iris$Species)
#' @export


pRank <- function(x,xSample){

  stopifnot(length(xSample) == nrow(x))

  x <- cbind.data.frame(x,"xSample"=xSample)%>%
    dplyr::group_by(xSample)%>%
    dplyr::mutate_all(.funs=rank.ub.average)%>%
    as.data.frame()%>%
    dplyr::select(-xSample)
  return(x)
}
