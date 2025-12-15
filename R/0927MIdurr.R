#setwd("C:\\Users\\Deng\\Dropbox\\MI R codes")
#' Multiple Imputation through Direct Use of Regularized Regression(DURR) for univariate missing data pattern
#' 
#' @param data A data frame or a matrix containing the incomplete data. Missing values are coded as NA.
#' @param family Specify the type of the missing data. Can be "gaussian", or "binary", or "poisson".
#' @param method The method to be used in regularized regression. Can be "lasso", "elastic net","adaptive lasso" or "blasso".
#' @param m Number of multiple imputations. The default is m=5.
#' @return n.missing Number of missing values in the data.
#' @return col.missing Which column contains the missing values.
#' @return impute The imputations(column number equals to m). Each column corresponds to one imputation.
#' @export
#' @examples
#' data(datagaussian)
#' data(databinary)
#' data(datapoisson)
#' MIdurr(data=datagaussian,method="lasso",family="gaussian")
#' MIdurr(data=datagaussian,method="blasso",family="gaussian")
#' MIdurr(data=databinary,method="elastic net",family="binary")
#' MIdurr(data=databinary,method="adaptive lasso",family="binary")
#' MIdurr(data=datapoisson,method="elastic net",family="poisson")
#' MIdurr(data=datapoisson,method="adaptive lasso",family="poisson")



MIdurr=function(data,family=c("gaussian"),method=c("lasso"),m=5)
{
realdata=as.matrix(data)
n=nrow(realdata)

colflag=as.numeric(which(colSums(is.na(realdata))!=0))
x1=realdata[,colflag]
nmis=sum(is.na(x1))
newobs<-realdata[,-colflag]
inpute=matrix(NA,n,m)

#tmp.type=length(unique(x1[!is.na(x1)]))/length(x1[!is.na(x1)])
#if(tmp.type>0.05){
#  type.missing=c("gaussian")
#} else if(length(unique(x1[!is.na(x1)]))==2){
#  type.missing=c("binary")
#} else {
#  type.missing=c("poisson")
#}
#
#family=type.missing



for(mm in 1:m){
  
  bootdata=realdata[sample(1:n,n,replace=TRUE),]
  fulldata=bootdata  
  comdata=fulldata[complete.cases(fulldata),]
  
  res<-comdata[,colflag]
  obs<-comdata[,-colflag]
  res<-as.matrix(res)
  obs<-as.matrix(obs)
  nnn<-length(res)
  
#######################################################
#Gaussian##############################################
#######################################################

if(family=="gaussian"){
##lasso
if(method=="lasso"){
  
  #model<- glmnet(obs,res,standardize=FALSE)
  cvob1=cv.glmnet(obs,res)
  predd<-predict(cvob1,obs,s="lambda.min")
  sse<-sqrt(sum((res-predd)^2)/nnn)
  pred<-predict(cvob1,newobs,s="lambda.min")
  
}

##elastic
if(method=="elastic net"){
  
  #model<- glmnet(obs,res,standardize=FALSE,alpha=0.5)
  cvob1=cv.glmnet(obs,res,alpha=0.5)
  predd<-predict(cvob1,obs,s="lambda.min")
  sse<-sqrt(sum((res-predd)^2)/nnn)
  pred<-predict(cvob1,newobs,s="lambda.min")
  
}

##adaptive lasso
if(method=="adaptive lasso"){
  
  betals=coef(cv.glmnet(obs,res,alpha=0))
  cvob1=cv.glmnet(obs,res,alpha=1,penalty.factor=1/abs(betals))
  predd<-predict(cvob1,obs,s="lambda.min")
  sse<-sqrt(sum((res-predd)^2)/nnn)
  pred<-predict(cvob1,newobs,s="lambda.min")
}

## blasso
if(method=="blasso"){
  model<-blasso.vs(res,obs,sig2=1,tau=1,phi=0.5,thin=2,sig2prior=c(0.01,0.01),tauprior=c(0.01,0.01),iters=5000,beta=rep(0.1,ncol(realdata)))
  bbb=apply(model$beta[3000:5000,],2,median)
  predd<-obs%*%bbb
  sse<-sqrt(sum((res-predd)^2)/nnn)
  pred<-newobs%*%bbb
}

##inpute missing values
mu4<-0
inpute[,mm]<-x1
for(i in 1:n){
  if(is.na(inpute[i,mm])){
    mu4<-pred[i]
    inpute[i,mm]<-rnorm(1,mu4,sd=sse)
  }
}

}


#######################################################
#Binary##############################################
#######################################################
if(family=="binary"){
  ##lasso
  if(method=="lasso"){
    
    #model<- glmnet(obs,res,standardize=FALSE)
    cvob1=cv.glmnet(obs,res,family="binomial")
    pred<-predict(cvob1,newobs,s="lambda.min")
    
  }
  
  ##elastic
  if(method=="elastic net"){
    
    #model<- glmnet(obs,res,standardize=FALSE,alpha=0.5)
    cvob1=cv.glmnet(obs,res,family="binomial",alpha=0.5)
    pred<-predict(cvob1,newobs,s="lambda.min")
    
  }
  
  ##adaptive lasso
  if(method=="adaptive lasso"){
    
    #model<-adalasso(obs,res)
    
    betals=coef(cv.glmnet(obs,res,family="binomial",alpha=0))
    
    # betals as numeric vector
    betals_vec <- as.numeric(betals)
    
    # Skip the intercept (usually the first element)
    penalty.factor <- 1 / abs(betals_vec[-1])   # now length == ncol(obs1)
    
    # Replace any non-finite values (Inf, NaN) with 1 for safety
    penalty.factor[!is.finite(penalty.factor)] <- 1
    cvob1=cv.glmnet(obs,res,alpha=1,penalty.factor=penalty.factor,family="binomial")
    pred<-predict(cvob1,newobs,s="lambda.min")
  }
  ##inpute missing values
  mu4<-0
  inpute[,mm]<-x1
  for(i in 1:n){
    if(is.na(inpute[i,mm])){
      mu4<-exp(pred[i])/(1+exp(pred[i]))
      inpute[i,mm]<-rbinom(1,1,mu4)
    }
  }
  
}




#######################################################
#Poisson##############################################
#######################################################
if(family=="poisson"){
  ##lasso
  if(method=="lasso"){
    
    #model<- glmnet(obs,res,standardize=FALSE)
    cvob1=cv.glmnet(obs,res,family="poisson")
    pred<-predict(cvob1,newobs,s="lambda.min")
    
  }
  
  ##elastic
  if(method=="elastic net"){
    
    #model<- glmnet(obs,res,standardize=FALSE,alpha=0.5)
    cvob1=cv.glmnet(obs,res,family="poisson",alpha=0.5)
    pred<-predict(cvob1,newobs,s="lambda.min")
    
  }
  
  ##adaptive lasso
  if(method=="adaptive lasso"){
    
    #model<-adalasso(obs,res)
    
    betals=coef(cv.glmnet(obs,res,family="poisson",alpha=0))
    # betals as numeric vector
    betals_vec <- as.numeric(betals)
    
    # Skip the intercept (usually the first element)
    penalty.factor <- 1 / abs(betals_vec[-1])   # now length == ncol(obs1)
    
    # Replace any non-finite values (Inf, NaN) with 1 for safety
    penalty.factor[!is.finite(penalty.factor)] <- 1
    cvob1=cv.glmnet(obs,res,alpha=1,penalty.factor=penalty.factor,family="poisson")
    pred<-predict(cvob1,newobs,s="lambda.min")
  }
  ##inpute missing values
  mu4<-0
  inpute[,mm]<-x1
  for(i in 1:n){
    if(is.na(inpute[i,mm])){
      mu4<-exp(pred[i])
      inpute[i,mm]<-rpois(1,mu4)
    }
  }
  
}


}  #end of mm



relist=list(n.missing=nmis,col.missing=colflag,impute=inpute)


return(relist)


}#end of function