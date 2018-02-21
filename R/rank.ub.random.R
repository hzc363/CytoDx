#' Percentile rank transformation of a vector (ties.method = "random")
#'
#' A function that performs the rank transformation of a vector (ties.method = "random").
#'
#' @param x A numeric vector.
#' @return Returns the percentile rank of each element.
#' @examples
#' rank.ub.average(1:10)
#' @export
rank.ub.random = function(x){
  n = length(unique(x))
  N = length(x)
  if(n==1){x=rep(NA,N)}
  if(n>1){x =(rank(x,ties.method="random")-1)/(N-1)*100}
  return(x)
}
