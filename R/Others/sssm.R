#' Build the second stage statistical model (sssm)
#'
#' A function that builds the second stage statistical model (sssm).
#'
#' @param fsp First stage prediction. Must be the object returned by fssm.glmnet or fssm.gbm.
#' @param xNew.Truth Optional. The dependent variable in the testing data if available.
#' @param type Type of prediction. Type of prediction required. Type "link"
#'   gives the linear predictors for "binomial", "multinomial", "poisson" or
#'   "cox" models; for "gaussian" models it gives the fitted values. Type
#'   "response" gives the fitted probabilities for "binomial" or "multinomial",
#'   fitted mean for "poisson" and the fitted relative-risk for "cox"; for
#'   "gaussian" type "response" is equivalent to type "link".
#' @return Returns a list. xNew.Pred contains the predicted y for the test data at the sample level.
#'   train.Data contains the trainig data and the predicted y for the training
#'   data at the sample level.  model contains the second stage statistical model. family specifies the type of regression.
#' @examples
#' # prepare data
#' x=model.matrix(~. , trainData[,1:7])
#' y=(trainData$Label=="aml")
#' xSample=trainData$SampleNumber
#' xNew=model.matrix(~. , testData[,1:7])
#' xSampleNew=testData$SampleNumber

#' # build the first stage statistical modeling
#' fsp =fssm.glmnet(x=x,y=y,xSample=xSample,xNew=xNew,reg=FALSE,
#'                  xSampleNew=xSampleNew,family="binomial")
#' # build the second stage statistical model
#' ssp = sssm(fsp)
#' @importFrom dplyr %>% group_by summarise_all left_join
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#' @export
sssm = function(fsp,xNew.Truth=NULL,type="response"){

  SP.train=fsp$train.Data%>%group_by(sample)%>%summarise_all(mean)
  SP.new=fsp$xNew.Pred%>%group_by(sample)%>%summarise_all(mean)


  y=SP.train[,grep(".Truth",colnames(SP.train))]
  y=as.matrix(y)
  if(fsp$family=="cox"){colnames(y)=c("time","status")}
  x=SP.train[,grep(".Pred",colnames(SP.train))]
  x =as.matrix(cbind("constant"=1,x))
  fit2 = glmnet(x=x,y=y,family=fsp$family,lambda=0)

  train.Data = predict(fit2,x,type=type)
  SP.train[,grep(".Pred",colnames(SP.train))]=as.data.frame(train.Data)
  names(SP.train)=colnames(fsp$train.Data)


  xNew = SP.new[,grep(".Pred",colnames(SP.new))]%>%as.matrix()
  xNew = cbind("constant"=1,xNew)
  xNew.Pred = predict(fit2,xNew,type=type)
  SP.new[,grep(".Pred",colnames(SP.new))] = as.data.frame(xNew.Pred)
  #SP.new=as.data.frame(SP.new)
  names(SP.new)=colnames(fsp$xNew.Pred)


  if(!is.null(xNew.Truth)){
    xNew.Truth = as.matrix(xNew.Truth)
    if(is.null(colnames(xNew.Truth))){colnames(xNew.Truth)=paste0("y",1:ncol(xNew.Truth))}
    LB = cbind.data.frame(fsp$xNew.Pred$sample,xNew.Truth)
    colnames(LB)=c("sample", paste0(colnames(xNew.Truth),".Truth"))
    LB = group_by(LB,sample)%>%summarise_all(unique)
    SP.new = left_join(SP.new,LB, by="sample")
  }

  return(list("xNew.Pred"=as.data.frame(SP.new),
              "train.Data"=as.data.frame(SP.train),
              "model"=fit2,
              "family"=fsp$family))
}
