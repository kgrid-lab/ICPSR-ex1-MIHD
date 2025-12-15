#' Multiple Imputation through Direct Use of Regularized Regression(DURR) in the presence of General missing patterns
#' 
#' @param data A data frame or a matrix containing the incomplete data with general missing patterns. Missing values are coded as NA.
#' @param family Specify the type of variables that have missing values. Can be either a single string or a vector of string with length of the number of variables having missing values. Each string can be "gaussian", or "binary", or "poisson".
#' @param method The method to be used in regularized regression. Can be "lasso", "elastic net", and "adaptive lasso".
#' @param m Number of multiple imputations. The default is \code{m=5}.
#' @param included If \code{TRUE}, GMdurr will include all observed values among m complete datasets.
#' @param burn A scalar giving the number of burns for the convergence of iterations. The default is \code{20}.
#' @return imp If \code{included=FALSE}, a list of components with the generated multiple imputations. The number of components equals the number of variables in the data set that have missing values. Each part of the list is a nmis[j] by m matrix of imputed values for variable j.
#' @return imp.included If \code{included=TRUE}, a list of components with the generated multiple imputations. The number of components equals the number of variables in the data set that have missing values. Each part of the list is a nrow(data) by m matrix of imputed values and observed values for variable j.
#' @export
#' @examples
#' data(GMdata)
#' GMdurr(data=GMdata)

GMdurr=function(data,m=5,family=rep("gaussian",3),method=c("lasso"),included=FALSE,burn=20){
ar1<-function(nn,rho){
  dudu<-matrix(rep(NA,nn*nn),ncol=nn,nrow=nn)
  for(i in 1:nn){
    for(j in 1:nn){
      dudu[i,j]<-rho^abs(i-j)
    }
  }
  return(dudu)
}
#############################################
#######Function begin here###################
#############################################
impute.dy=function(W_mis,theta,family=c("gaussian")){
  W_mis=as.matrix(W_mis)
  if(family=="gaussian"){
    z_mis=rnorm(n=nrow(W_mis),mean=as.matrix(cbind(1,W_mis))%*%theta[-length(theta)],sd=sqrt(theta[length(theta)]))   
  }
  if(family=="binary"){
    nu=exp(as.matrix(cbind(1,W_mis))%*%theta)
    z_mis=rbinom(n=nrow(W_mis),size=1,prob=nu/(1+nu))
  }
  if(family=="poisson"){
    nu=exp(as.matrix(cbind(1,W_mis))%*%theta)
    z_mis=rpois(n=nrow(W_mis),lambda=nu)
  }
  return(z_mis)
}

act.set.dy=function(z_obs,W_obs,family=c("gaussian"),method=c("lasso")){
  nochange=FALSE
  p=length(z_obs)
  if(family=="gaussian"){
    if(method=="lasso"){
      model<- glmnet(W_obs,z_obs,pmax=p-1,standardize=FALSE)
      cvres<-cv.glmnet(W_obs,z_obs,pmax=p-1,standardize=FALSE)
      fits<-coef(model, s=cvres$lambda.min)
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
    if(method=="elastic_net"){
      model<- glmnet(W_obs,z_obs,pmax=p-1,standardize=FALSE,alpha=0.5)
      cvres<-cv.glmnet(W_obs,z_obs,pmax=p-1,standardize=FALSE,alpha=0.5)
      fits<-coef(model, s=cvres$lambda.min)
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
    if(method=="adaptive_lasso"){
      model<- adalasso(W_obs,z_obs)
      fits<-model$coefficients.adalasso
      aset=try(as.vector(which(fits!=0)),silent=TRUE)
    }
  }
  
  
  
  if(family=="binary"){
    if(method=="lasso"){
      model<- glmnet(W_obs,z_obs,family="binomial",pmax=p-1,standardize=FALSE)
      cvres<-cv.glmnet(W_obs,z_obs,family="binomial",pmax=p-1,standardize=FALSE)
      fits<-coef(model, s=cvres$lambda.min)
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
    if(method=="elastic_net"){
      model<- glmnet(W_obs,z_obs,family="binomial",pmax=p-1,standardize=FALSE,alpha=0.5)
      cvres<-cv.glmnet(W_obs,z_obs,family="binomial",pmax=p-1,standardize=FALSE,alpha=0.5)
      fits<-coef(model, s=cvres$lambda.min)
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
    if(method=="adaptive_lasso"){
      theridge.cv<-cv.glmnet(W_obs,z_obs,family="binomial",pmax=p-1,standardize=FALSE,alpha=0) ## first stage ridge
      
      ## Second stage weights from the coefficients of the first stage
      bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
      if(all(bhat==0)){
        ## if bhat is all zero then assign very close to zero weight to all.
        ## Amounts to penalizing all of the second stage to zero.
        bhat<-rep(.Machine$double.eps*2,length(bhat))
      }
      adpen<-(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
      
      ## Second stage lasso (the adaptive lasso)
      thelasso.cv<-cv.glmnet(W_obs,z_obs,family="binomial",pmax=p-1,standardize=FALSE,alpha=1,
                             exclude=which(bhat==0),
                             penalty.factor=adpen)
      ## Extract resulting coefs
      fits<-coef(thelasso.cv,s="lambda.min")
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
  }
  
  
  if(family=="poisson"){
    if(method=="lasso"){
      model<- glmnet(W_obs,z_obs,family="poisson",pmax=p-1,standardize=FALSE)
      cvres<-cv.glmnet(W_obs,z_obs,family="poisson",pmax=p-1,standardize=FALSE)
      fits<-coef(model, s=cvres$lambda.min)
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
    if(method=="elastic_net"){
      model<- glmnet(W_obs,z_obs,family="poisson",pmax=p-1,standardize=FALSE,alpha=0.5)
      cvres<-cv.glmnet(W_obs,z_obs,family="poisson",pmax=p-1,standardize=FALSE,alpha=0.5)
      fits<-coef(model, s=cvres$lambda.min)
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
    if(method=="adaptive_lasso"){
      theridge.cv<-cv.glmnet(W_obs,z_obs,family="poisson",pmax=p-1,standardize=FALSE,alpha=0) ## first stage ridge
      
      ## Second stage weights from the coefficients of the first stage
      bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
      if(all(bhat==0)){
        ## if bhat is all zero then assign very close to zero weight to all.
        ## Amounts to penalizing all of the second stage to zero.
        bhat<-rep(.Machine$double.eps*2,length(bhat))
      }
      adpen<-(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
      
      ## Second stage lasso (the adaptive lasso)
      thelasso.cv<-cv.glmnet(W_obs,z_obs,family="poisson",pmax=p-1,standardize=FALSE,alpha=1,
                             exclude=which(bhat==0),
                             penalty.factor=adpen)
      ## Extract resulting coefs
      fits<-coef(thelasso.cv,s="lambda.min")
      aset=as.matrix(summary(fits))[,1]-1
      if(aset[1]==0){aset=aset[-1]}  
    }
  }   
  if(length(aset)==0){nochange=TRUE}
  return(list(aset=aset,nochange=nochange))
}

durr.pred=function(z_obs,W_obs,W_mis,family=c("gaussian"),method=c("lasso")){
  nochange=FALSE
  p=length(z_obs)
  if(family=="gaussian"){
    if(method=="lasso"){
      cvob1=cv.glmnet(W_obs,z_obs)
      predd<-predict(cvob1,W_obs,s="lambda.min")
      sse<-sqrt(sum((z_obs-predd)^2)/p)
      pred<-predict(cvob1,W_mis,s="lambda.min")
    }
    if(method=="elastic_net"){
      cvob1=cv.glmnet(W_obs,z_obs,alpha=0.5)
      predd<-predict(cvob1,W_obs,s="lambda.min")
      sse<-sqrt(sum((z_obs-predd)^2)/p)
      pred<-predict(cvob1,W_mis,s="lambda.min")
    }
    if(method=="adaptive_lasso"){
      theridge.cv<-cv.glmnet(W_obs,z_obs,standardize=FALSE,alpha=0) ## first stage ridge
      
      ## Second stage weights from the coefficients of the first stage
      bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
      if(all(bhat==0)){
        ## if bhat is all zero then assign very close to zero weight to all.
        ## Amounts to penalizing all of the second stage to zero.
        bhat<-rep(.Machine$double.eps*2,length(bhat))
      }
      adpen<-(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
      
      ## Second stage lasso (the adaptive lasso)
      thelasso.cv<-cv.glmnet(W_obs,z_obs,standardize=FALSE,alpha=1,
                             exclude=which(bhat==0),
                             penalty.factor=adpen)
      predd<-predict(thelasso.cv,W_obs,s="lambda.min")
      sse<-sqrt(sum((z_obs-predd)^2)/p)
      pred<-predict(thelasso.cv,W_mis,s="lambda.min")
    }
    impute=rnorm(n=length(pred),mean=pred,sd=sse)
  } 
  
  if(family=="binary"){
    if(method=="lasso"){
      cvob1=cv.glmnet(W_obs,z_obs,family="binomial")
      pred<-predict(cvob1,W_mis,s="lambda.min")
    }
    if(method=="elastic_net"){
      cvob1=cv.glmnet(W_obs,z_obs,family="binomial",alpha=0.5)
      pred<-predict(cvob1,W_mis,s="lambda.min")
    }
    if(method=="adaptive_lasso"){
      theridge.cv<-cv.glmnet(W_obs,z_obs,family="binomial",standardize=FALSE,alpha=0) ## first stage ridge
      
      ## Second stage weights from the coefficients of the first stage
      bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
      if(all(bhat==0)){
        ## if bhat is all zero then assign very close to zero weight to all.
        ## Amounts to penalizing all of the second stage to zero.
        bhat<-rep(.Machine$double.eps*2,length(bhat))
      }
      adpen<-(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
      
      ## Second stage lasso (the adaptive lasso)
      thelasso.cv<-cv.glmnet(W_obs,z_obs,family="binomial",standardize=FALSE,alpha=1,
                             exclude=which(bhat==0),
                             penalty.factor=adpen)
      pred<-predict(thelasso.cv,W_mis,s="lambda.min")
    }
    impute=rbinom(n=length(pred),size=1,prob=exp(pred)/(1+exp(pred)))
  }
  
  
  if(family=="poisson"){
    if(method=="lasso"){
      cvob1=cv.glmnet(W_obs,z_obs,family="poisson")
      pred<-predict(cvob1,W_mis,s="lambda.min")
    }
    if(method=="elastic_net"){
      cvob1=cv.glmnet(W_obs,z_obs,family="poisson",alpha=0.5)
      pred<-predict(cvob1,W_mis,s="lambda.min")
    }
    if(method=="adaptive_lasso"){
      theridge.cv<-cv.glmnet(W_obs,z_obs,family="poisson",standardize=FALSE,alpha=0) ## first stage ridge
      
      ## Second stage weights from the coefficients of the first stage
      bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
      if(all(bhat==0)){
        ## if bhat is all zero then assign very close to zero weight to all.
        ## Amounts to penalizing all of the second stage to zero.
        bhat<-rep(.Machine$double.eps*2,length(bhat))
      }
      adpen<-(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
      
      ## Second stage lasso (the adaptive lasso)
      thelasso.cv<-cv.glmnet(W_obs,z_obs,family="poisson",standardize=FALSE,alpha=1,
                             exclude=which(bhat==0),
                             penalty.factor=adpen)
      pred<-predict(thelasso.cv,W_mis,s="lambda.min")
    }
    impute=rpois(n=length(pred),lambda=exp(pred))
  }
  
  return(list(impute=impute,nochange=nochange))
}




  data=as.matrix(data)
  m.final=m
  m=m+burn
  n=nrow(data)
  p=ncol(data)
  mis.pat=is.na(data) #missing pattern
  nmis=apply(mis.pat,2,sum) #An array containing the number of missing observations per column
  ind.col=as.vector(which(nmis!=0)) #index of which column has missingness  
  colna=colnames(data)[ind.col] #names of columns that have missingness
  l=length(ind.col)
  Z0=data
  ini.mean=apply(data,2,function(o) mean(o,na.rm=TRUE))
  ind.row=list() #a list of l vectors of which rows are missing for variable j
  imp=list() #a list imputations consider the initial mean value to be the first imputation (will delete in the future)
  
  for(ll in 1:l){ind.row[[colna[ll]]]=as.vector(which(is.na(data[,ind.col[ll]])==TRUE))
  imp[[colna[ll]]]=matrix(rep(ini.mean[ind.col[ll]],length(which(is.na(data[,ind.col[ll]])==TRUE))),ncol=1)
  rownames(imp[[colna[ll]]])=as.list(ind.row[[colna[ll]]])
  Z0[ind.row[[colna[ll]]],ind.col[ll]]=ini.mean[ind.col[ll]]
  } 
  
  for(mm in 1:m){
    for(j in 1:l){
      W=Z0[,-ind.col[j]]
      boos=sample(1:n,size=n,replace=TRUE)
      for(k in 1:length(ind.row[[j]])){
        boos=boos[boos!=ind.row[[j]][k]]  #observed row in after bootstrap 
      }
      ast=boos
      
      impute=NULL  
      #impute=tryCatch(durr.pred(z_obs=data[ast,ind.col[j]],W_obs=W[ast,],W_mis=W[-ind.row[[j]],],family=family[j],method=method)$impute,function(e) NULL)
      try(impute<-durr.pred(z_obs=Z0[ast,ind.col[j]],W_obs=W[ast,],W_mis=W[ind.row[[j]],],family=family[j],method=method)$impute,silent=TRUE)
      impute=unname(impute)
      if(length(impute)==0){
        imp[[j]]=cbind(imp[[j]],imp[[j]][,ncol(imp[[j]])])
        next}
      
      imp[[j]]=cbind(imp[[j]],impute)
      Z0[ind.row[[j]],ind.col[[j]]]=imp[[j]][,mm+1]
    }
  }
  
  for (name in names(imp)){imp[[name]]=imp[[name]][,-(1:(burn+1))]}
  if(included==FALSE){
    return(list(imp=imp))  
  } else{
    imp.included=list()
    for(ll in 1:l){
      imp.included[[colna[ll]]]=matrix(rep(data[,ind.col[ll]],m.final),ncol=m.final)
      imp.included[[colna[ll]]][ind.row[[colna[ll]]],]=imp[[colna[ll]]]
    }
    return(list(imp.included=imp.included))   
  }  
}

