meanUnique = function(x){
  if(is.numeric(x)){
    return(mean(x,na.rm = T))
  }else(return(unique(x)))
}

tssm.glm.matrix= function(x,y,xSample,family,type1="response",type2="response",
                          parallelCore=1,reg=FALSE,newoffset=NULL,...){
  y = as.matrix(y)
  if(parallelCore>1){doParallel::registerDoParallel(parallelCore)}


  if(reg==TRUE){
    fit <- glmnet::cv.glmnet(x=x,y=y,family=family,parallel=(parallelCore>1),...)
  }else{
    fit <- glmnet::glmnet(x=x,y=y,family=family,lambda = 0,...)
  }

  preds = predict(fit,newx = x,type=type1,newoffset=newoffset)
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
  names(SP.train)=colnames(train.Data)

  return(list("train.Data1"=train.Data,"model1"=fit,"train.Data2"=as.data.frame(SP.train),
              "model2"=fit2,"family"=family,"method"="glm","type1"=type1,"type2"=type2))
}



tssm.pred.matrix = function(fit,xNew,xSampleNew){

  if(fit$method%in%c("glm","nnet")){
    preds = predict(fit$model1,xNew,type=fit$type1)
  }else if(fit$method%in%c("gbm")){
    preds = predict(fit$model1,xNew,type=fit$type1,n.trees =fit$model1$n.trees)
    preds = as.matrix(preds)
  }else{print("Please use output from tssm.glm or tssm.gmb.")}

  xNew.Pred = cbind.data.frame("sample"=xSampleNew,preds)
  colnames(xNew.Pred)=c("sample",paste0("y",".Pred.",colnames(preds)))
  SP.new=xNew.Pred%>%group_by(sample)%>%summarise_all(mean)

  xNew = SP.new[,grep(".Pred",colnames(SP.new))]%>%as.matrix()
  xNew = cbind("constant"=1,xNew)
  xNew.Pred2 = predict(fit$model2,xNew,type=fit$type2)
  SP.new[,grep(".Pred",colnames(SP.new))] = as.data.frame(xNew.Pred2)
  #SP.new=as.data.frame(SP.new)
  names(SP.new)=colnames(xNew.Pred)

  return(list("xNew.Pred1"=xNew.Pred,
              "xNew.Pred2"=as.data.frame(SP.new)))

}
