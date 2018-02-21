#' Build the CytoDx model
#'
#' A function that builds the CytoDx model.
#'
#' @param x The marker profile of cells pooled from all samples. Each row is a cell, each column is
#'   a marker.
#' @param y The clinical outcomes associated with samples to which cells belong. Length must be equal to nrow(x). For
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
#' @param xSample A vector specifying which sample each cell belongs to. Length
#'   must equal to nrow(x).
#' @param family Response type. Must be one of the following:
#'   "gaussian","binomial","poisson","multinomial","cox","mgaussian"
#' @param type1 Type of first level prediction. Type of prediction required.
#'   Type "link" gives the linear predictors for "binomial", "multinomial",
#'   "poisson" or "cox" models; for "gaussian" models it gives the fitted
#'   values. Type "response" gives the fitted probabilities for "binomial" or
#'   "multinomial", fitted mean for "poisson" and the fitted relative-risk for
#'   "cox"; for "gaussian" type "response" is equivalent to type "link".
#' @param type2 Type of second level prediction.
#' @param parallelCore The number of core to be used. Only used when reg is
#'   TRUE.
#' @param reg If elestic net regularization will be used.
#' @param ... Other parameters to be passed into the glmnet or the cv.glmnet
#'   function in the glmnet package.
#' @return Returns a list. train.Data.cell contains the trainig data and the
#'   predicted y for the training data at the cell level. model.cell contains the
#'   cell stage statistical model. Data.sample contains the trainig data and the
#'   predicted y for the training data at the sample level. model.sample contains the
#'   sample stage statistical model. family specifies the regression type.
#'   method specifies the type of learning method. type.cell is the type of cell
#'   level prediction. type.sample is the type of sample level prediction.
#' @examples
#' library(CytoDx)
#' # Find the table containing fcs file names in CytoDx package
#' path=system.file("extdata",package="CytoDx")
#' # read the table
#' fcs_info = read.csv(file.path(path,"fcs_info.csv"))
#' # Specify the path to the cytometry files
#' fn = file.path(path,fcs_info$fcsName)
#' # Read cytometry files using fcs2DF function
#' train_data = fcs2DF(fcsFiles=fn,
#'                     y=fcs_info$Label,
#'                     assay="FCM",
#'                     b=1/150,
#'                     excludeTransformParameters=
#'                       c("FSC-A","FSC-W","FSC-H","Time"))
#' # build the model
#' fit =CytoDx.fit(x=as.matrix(train_data[,1:7]),
#'                 y=train_data$y,
#'                 xSample = train_data$xSample,
#'                 reg=FALSE,
#'                 family="binomial")
#' # check accuracy for training data
#' pred = CytoDx.pred(fit,
#'                    xNew=as.matrix(train_data[,1:7]),
#'                    xSampleNew=train_data$xSample)
#'
#' boxplot(pred$xNew.Pred.sample$y.Pred.s0~
#'           fcs_info$Label)
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats predict
#' @importFrom dplyr %>% group_by summarise_all left_join
#' @export

CytoDx.fit= function(x,y,xSample,family,type1="response",type2="response",
                          parallelCore=1,reg=FALSE,...){
  y = as.matrix(y)
  xSample=as.character(xSample)
  if(parallelCore>1){doParallel::registerDoParallel(parallelCore)}


  if(reg==TRUE){
    fit <- glmnet::cv.glmnet(x=x,y=y,family=family,
                             parallel=(parallelCore>1),...)
  }else{
    fit <- glmnet::glmnet(x=x,y=y,family=family,lambda = 0,...)
  }

  preds = predict(fit,newx = x,type=type1)
  train.Data = cbind.data.frame(xSample,y,preds)
  y=as.matrix(y)
  colnames(train.Data)=c("sample",paste0("y",1:ncol(y),".Truth"),
                         paste0("y",".Pred.",colnames(preds)))

  SP.train=train.Data%>%group_by(sample)%>%summarise_all(meanUnique)


  y=SP.train[,grep(".Truth",colnames(SP.train))]
  y=as.matrix(y)
  if(family=="cox"){colnames(y)=c("time","status")}
  x=SP.train[,grep(".Pred",colnames(SP.train))]
  x =as.matrix(cbind("constant"=1,x))
  fit2 = glmnet::glmnet(x=x,y=y,family=family,lambda=0)

  train.Data2 = predict(fit2,x,type=type2)
  SP.train[,grep(".Pred",colnames(SP.train))]=as.data.frame(train.Data2)
  SP.train=as.data.frame(SP.train)
  colnames(SP.train)=colnames(train.Data)
  row.names(SP.train)=SP.train$sample
  SP.train = SP.train[unique(xSample),]
  return(list("train.Data.cell"=train.Data,"model.cell"=fit,"train.Data.sample"=SP.train,
              "model.sample"=fit2,"family"=family,"type.cell"=type1,"type.sample"=type2))
}
