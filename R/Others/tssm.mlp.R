#' Build the two stage statistical model (tssm) using neural network
#'
#' A function that builds the two stage statistical model (tssm) using neural network.
#'
#' @param x A matrix of training data. Each row is a record, each column is a
#'   independent variable included in the model. The records in all samples
#'   should be pasted into a single matrix.
#' @param y The dependent variable in the training data. The records that
#'   belongs to the same sample should have the same dependent variable. For
#'   family="binomial" should be either a factor with two levels, or a
#'   two-column matrix of counts or proportions (the second column is treated as
#'   the target class; for a factor, the last level in alphabetical order is the
#'   target class). For family="multinomial", can be a nc>=2 level factor, or a
#'   matrix with nc columns of counts or proportions. For either "binomial" or
#'   "multinomial", if y is presented as a vector, it will be coerced into a
#'   factor. For family="cox", y should be a two-column matrix with columns
#'   named 'time' and 'status'. The latter is a binary variable, with '1'
#'   indicating death, and '0' indicating right censored. The function Surv() in
#'   package survival produces such a matrix. For family="mgaussian", y is a
#'   matrix of quantitative responses.
#' @param xSample A vector of sample ID in the training data. All the records
#'   that belongs to the same sample should have the same sample ID.
#' @param family Response type. Must be one of the following:
#'   "gaussian","binomial".
#' @param type1 Type of first level prediction. Type of prediction required.
#'   Type "link" gives the linear predictors for "binomial", "multinomial",
#'   "poisson" or "cox" models; for "gaussian" models it gives the fitted
#'   values. Type "response" gives the fitted probabilities for "binomial" or
#'   "multinomial", fitted mean for "poisson" and the fitted relative-risk for
#'   "cox"; for "gaussian" type "response" is equivalent to type "link".
#' @param type2 Type of second level prediction.
#' @param ... Other parameters to be passed into the nnet function in the nnet library.
#' @return Returns a list. train.Data1 contains the trainig data and the
#'   predicted y for the training data at the first level. model1 contains the
#'   first stage statistical model. Data2 contains the trainig data and the
#'   predicted y for the training data at the second level. model2 contains the
#'   second stage statistical model. family specifies the regression type.
#'   method specifies the type of learning method. type1 is the type of first
#'   level prediction. type2 is the type of second level prediction.
#' @examples
#' # prepare data
#' x=model.matrix(~. , trainData[,1:7])
#' y=(trainData$Label=="aml")
#' xSample=trainData$SampleNumber

#' # build the tssm
#' fit =tssm.mlp(x=x,y=y,xSample=xSample,reg=FALSE,family="binomial",linout=T)
#' @importFrom RSNNS mlp
#' @importFrom stats predict
#' @export
tssm.mlp = function(x,y,xSample,family,type1="raw",type2="response",size = 20,linout=TRUE,...){


  fit <- mlp(x,y,size=size,linOut=linout,...)

  preds = predict(fit,newx = x,type=type1)
  train.Data = cbind.data.frame(xSample,y,preds)
  y=as.matrix(y)
  colnames(train.Data)=c("sample",paste0("y",1:ncol(y),".Truth"),paste0("y",1:ncol(preds),".Pred"))

  SP.train=train.Data%>%group_by(sample)%>%summarise_all(mean)

  y=SP.train[,grep(".Truth",colnames(SP.train))]
  y=as.matrix(y)
  if(family=="cox"){colnames(y)=c("time","status")}
  x=SP.train[,grep(".Pred",colnames(SP.train))]
  x =as.matrix(cbind("constant"=1,x))
  fit2 = glmnet::glmnet(x=x,y=y,family=family,lambda=0)

  train.Data2 = predict(fit2,x,type=type2)
  SP.train[,grep(".Pred",colnames(SP.train))]=as.data.frame(train.Data2)
  names(SP.train)=colnames(train.Data)

  return(list("train.Data1"=train.Data,"model1"=fit,"train.Data2"=as.data.frame(SP.train),
              "model2"=fit2,"family"=family,"method"="nnet","type1"=type1,"type2"=type2))
}

