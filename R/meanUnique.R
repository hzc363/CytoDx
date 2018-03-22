#' Calulate mean or take unique elements of a vector
#'
#' A function that calulate mean or take unique elements of a vector.
#'
#' @param x a vector
#' @return If x is numeric, returns the mean. Otherwise, returns the unique elements of x.
#' @examples
#' x = 1:5
#' meanUnique(x)
#' x=c("a","a","b")
#' meanUnique(x)
#' @export

meanUnique = function(x){
  if(is.numeric(x)){
    return(mean(x,na.rm = TRUE))
  }else(return(unique(x)))
}
