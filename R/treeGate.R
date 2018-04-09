#' Use decision tree to find a group of cells that are associated with clinical
#' outcome.
#'
#' A function that sse decision tree to find a group of cells that are
#' associated with clinical outcome.
#'
#' @param P The predicted association of each cell with a clinical outcome.
#' @param x The marker profile of each cell. Each row is a cell, each column is
#'   a marker. Must have length(P) rows.
#' @param ... Other parameters to be passed into the rpart function
#' @return Returns a object created by rpart function. Also plots a graph of
#'   decision tree.
#' @examples
#' # Find the table containing fcs file names in CytoDx package
#' path=system.file("extdata",package="CytoDx")
#' # read the table
#' fcs_info <- read.csv(file.path(path,"fcs_info.csv"))
#' # Specify the path to the cytometry files
#' fn <- file.path(path,fcs_info$fcsName)
#' # Read cytometry files using fcs2DF function
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
#' # Find the associated population using treeGate
#' TG <- treeGate(P = fit$train.Data.cell$y.Pred.s0,
#'               x= train_data[,1:7])
#' @importFrom rpart rpart rpart.control
#' @importFrom rpart.plot rpart.plot
#' @export
treeGate <- function(P,x,...){
  #if(is.null(names(x))){names(x)=paste0("sample",1:length(x))}
  #xN = sapply(x,nrow)
  #x = Reduce(rbind,x)

  colnames(x) <- gsub("[^[:alnum:]]", "_",colnames(x))
  #P = fit$train.Data1$y.Pred.s0
  x <- cbind.data.frame("P"=P,x)
  fit <- rpart(P~.,method="anova", data=x,...)
  rpart.plot(fit)
  return(fit)
}

