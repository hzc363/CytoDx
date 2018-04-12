#' Make prediction using the CytoDx model
#'
#' A function that makes prediction using the CytoDx model.
#'
#' @param fit The two stage statistical model. Must be the object returned by CytoDx.fit.
#' @param xNew The marker profile of cells pooled from all new samples. Each row is a cell, each column is
#'   a marker.
#' @param xSampleNew A vector specifying which sample each cell belongs to. Length
#'   must equal to nrow(xNew).
#' @return Returns a list. xNew.Pred1 contains the predicted y for the new data at the cell level.
#'   xNew.Pred2 contains the predicted y for the new data at the sample level.
#' @examples
#' # Find the table containing fcs file names in CytoDx package
#' path <- system.file("extdata",package="CytoDx")
#' # read the table
#' fcs_info <- read.csv(file.path(path,"fcs_info.csv"))
#' # Specify the path to the cytometry files
#' fn <- file.path(path,fcs_info$fcsName)
# Read cytometry files using fcs2DF function
#' train_data <- fcs2DF(fcsFiles=fn,
#'                     y=fcs_info$Label,
#'                     assay="FCM",
#'                     b=1/150,
#'                     excludeTransformParameters=
#'                       c("FSC-A","FSC-W","FSC-H","Time"))
#' # build the model
#' fit <- CytoDx.fit(x=as.matrix(train_data[,1:7]),
#'                 y=train_data$y,
#'                 xSample = train_data$xSample,
#'                 reg=FALSE,
#'                 family="binomial")
#' # check accuracy for training data
#' pred <- CytoDx.pred(fit,
#'                    xNew=as.matrix(train_data[,1:7]),
#'                    xSampleNew=train_data$xSample)
#'
#' boxplot(pred$xNew.Pred.sample$y.Pred.s0~
#'           fcs_info$Label)
#'
#' @importFrom dplyr %>% group_by summarise_all left_join
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#' @export

CytoDx.pred <- function(fit,xNew,xSampleNew){

  stopifnot(length(xSampleNew) == nrow(xNew))

  xSampleNew <- as.character(xSampleNew)

  preds <- predict(fit$model.cell,xNew,type=fit$type.cell)


  xNew.Pred <- cbind.data.frame("sample"=xSampleNew,preds)
  colnames(xNew.Pred) <- c("sample",paste0("y",".Pred.",colnames(preds)))
  SP.new <- xNew.Pred%>%group_by(sample)%>%summarise_all(mean)

  xNew <- SP.new[,grep(".Pred",colnames(SP.new))]%>%as.matrix()
  xNew <- cbind("constant"=1,xNew)
  xNew.Pred2 <- predict(fit$model.sample,xNew,type=fit$type.sample)
  SP.new[,grep(".Pred",colnames(SP.new))] <- as.data.frame(xNew.Pred2)

  SP.new <- as.data.frame(SP.new)
  colnames(SP.new) <- colnames(xNew.Pred)
  row.names(SP.new) <- SP.new$sample

  SP.new  <-  SP.new[unique(xSampleNew),]
  return(list("xNew.Pred.cell"=xNew.Pred,
              "xNew.Pred.sample"=SP.new))

}

