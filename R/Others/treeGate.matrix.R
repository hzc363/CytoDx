#' Use decision tree to find a group of records that are predictive.
#'
#' A function that use decision tree to find a group of records that are predictive.
#'
#' @param x Must be the same list as the x argument in tssm.glm.
#' @param pred The predicted association between each cell and the clinical outcome.
#' @param ... Other parameters to be passed into the rpart function
#' @return Returns a object created by rpart function. Also plots a graph of decision tree.
#' @importFrom rpart rpart rpart.control
#' @importFrom rpart.plot rpart.plot
#' @export
#'
treeGate.matrix = function(x,pred,...){
  colnames(x)=gsub("[^[:alnum:]]", "_",colnames(x))
  P = pred #fit$train.Data1$y1.Pred
  x = cbind.data.frame("P"=P,x)
  fit <- rpart(P~.,method="anova", data=x,...)
  rpart.plot(fit)
  return(fit)
}
