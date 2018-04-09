#' Percentile rank transformation of a vector
#'
#' A function that performs the Percentile rank transformation of a vector
#'
#' @param x A numeric vector.
#' @return Returns the percentile rank of each element.
#' @examples
#' rank.ub.average(1:10)
#' @export
rank.ub.average <- function(x){
  n <- length(unique(x))
  N <- length(x)
  if(n==1){x <- rep(NA,N)}
  if(n>1){x  <- (rank(x,ties.method="average")-1)/(N-1)*100}
  return(x)
}
