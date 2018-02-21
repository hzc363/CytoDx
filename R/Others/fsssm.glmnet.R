#' Build the first stage statistical model (fssm) using regularized regression
#'
#' A function that builds the first stage statistical model (fssm) using
#' regularized regression.
#' @param x A matrix of training data. Each row is a record, each column is a
#'   independent variable included in the model. The records in all samples
#'   should be pasted into a single matrix.
#' @param y The dependent variable in the training data. The records
#'   that belongs to the same sample should have the same dependent variable.
#'   For family="binomial" should be either a factor with two levels, or a
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
#' @param xNew Optional. A matrix of testing data. Each row is a record, each
#'   column is a independent variable included in the model. The records in all
#'   samples should be pasted into a single matrix.
#' @param xSampleNew A vector of sample ID in the testing data. All the records
#'   that belongs to the same sample should have the same sample ID.
#' @param family Response type. Must be one of the following:
#'   "gaussian","binomial","poisson","multinomial","cox","mgaussian"
#' @param type Type of prediction. Type of prediction required. Type "link"
#'   gives the linear predictors for "binomial", "multinomial", "poisson" or
#'   "cox" models; for "gaussian" models it gives the fitted values. Type
#'   "response" gives the fitted probabilities for "binomial" or "multinomial",
#'   fitted mean for "poisson" and the fitted relative-risk for "cox"; for
#'   "gaussian" type "response" is equivalent to type "link".
#' @param parallelCore The number of core to be used. Only used when reg is
#'   TRUE.
#' @param model The first stage glmnet model to be used. If NULL, a model will
#'   be built using training data.
#' @param reg If elestic net regularization will be used.
#' @param ... Other parameters to be passed into the glmnet or the cv.glmnet
#'   function in the glmnet package.
#' @return Returns a list. xNew.Pred contains the predicted y for the test data.
#'   train.Data contains the trainig data and the predicted y for the training
#'   data.  model contains the first stage statistical model. family specifies
#'   the regression type.
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
#' @importFrom doParallel registerDoParallel
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats predict
#' @export
fssm.glmnet = function(x,y,xSample,xNew=NULL,xSampleNew=NULL,family,type="response",
                       parallelCore=1,model=NULL,reg=FALSE,...){

  if(parallelCore>1){doParallel::registerDoParallel(parallelCore)}

  if(is.null(model)){
    if(reg==TRUE){
      fit <- glmnet::cv.glmnet(x=x,y=y,family=family,parallel=(parallelCore>1),...)
    }else{
      fit <- glmnet::glmnet(x=x,y=y,family=family,lambda = 0,...)
    }
  }else{fit=model}


  # put preds into the right scale
  if(is.null(xNew)){xNew=x;xSampleNew=xSample}
  preds = predict(fit,newx = x,type=type)
  train.Data = cbind.data.frame(xSample,y,preds)
  y=as.matrix(y)
  colnames(train.Data)=c("sample",paste0("y",1:ncol(y),".Truth"),paste0("y",1:ncol(preds),".Pred"))

  preds = predict(fit,xNew,type=type)
  xNew.Pred = cbind.data.frame("sample"=xSampleNew,preds)
  colnames(xNew.Pred)=c("sample",paste0("y",1:ncol(preds),".Pred"))

  return(list("xNew.Pred"=xNew.Pred,"train.Data"=train.Data,"model"=fit,"family"=family))
}

