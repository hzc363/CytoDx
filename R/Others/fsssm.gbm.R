

#' Build the first stage statistical model (fssm) using gradient boosting
#' machine.
#'
#' A function that builds the first stage statistical model (fssm) using
#' gradient boosting machine..
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
#' @param xNew Optional. A matrix of testing data. Each row is a record, each
#'   column is a independent variable included in the model. The records in all
#'   samples should be pasted into a single matrix.
#' @param xSampleNew A vector of sample ID in the testing data. All the records
#'   that belongs to the same sample should have the same sample ID.
#' @param family Response type. Must be one of the following:
#'   "gaussian","binomial","poisson","multinomial","cox","mgaussian"
#' @param type Type of prediction. If "link", the predictions are on the scale
#'   of f(x). For example, for the Bernoulli loss the returned value is on the
#'   log odds scale, poisson loss on the log scale, and coxph is on the log
#'   hazard scale.If type="response" then gbm converts back to the same scale as
#'   the outcome. Currently the only effect this will have is returning
#'   probabilities for bernoulli and expected counts for poisson. For the other
#'   distributions "response" and "link" return the same.
#' @param model The first stage glmnet model to be used. If NULL, a model will
#'   be built using training data.
#' @param ... Other parameters to be passed into the gbm.fit
#'   function in the gbm package.
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
#' fsp =fssm.gbm(x=x,y=y,xSample=xSample,xNew=xNew,
#'                  xSampleNew=xSampleNew,family="binomial")
#' @importFrom gbm gbm.fit
#' @importFrom stats predict
#' @export
fssm.gbm = function(x,y,xSample,xNew=NULL,xSampleNew=NULL,family,
                    type="response",model=NULL,...){
  familyConv=c("gaussian" , "bernoulli" ,"poisson","multinomial","coxph" )
  names(familyConv) = c("gaussian", "binomial", "poisson", "multinomial", "cox" )

  if(is.null(model)){
    fit <- gbm::gbm.fit(x,y,distribution =familyConv[family],...)
  }else{fit=model}

  # put preds into the right scale
  y=as.matrix(y)
  if(is.null(xNew)){xNew=x;xSampleNew=xSample}
  preds = predict(fit,newx = x,type=type,n.trees = fit$n.trees)
  train.Data = cbind.data.frame(xSample,y,preds)
  preds=as.matrix(preds)
  colnames(train.Data)=c("sample",paste0("y",1:ncol(y),".Truth"),paste0("y",1:ncol(preds),".Pred"))

  preds = predict(fit,xNew,type=type,n.trees = fit$n.trees)
  xNew.Pred = cbind.data.frame("sample"=xSampleNew,preds)
  preds=as.matrix(preds)
  colnames(xNew.Pred)=c("sample",paste0("y",1:ncol(preds),".Pred"))

  return(list("xNew.Pred"=xNew.Pred,"train.Data"=train.Data,"model"=fit,"family"=family))
}

