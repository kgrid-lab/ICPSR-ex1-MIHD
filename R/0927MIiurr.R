#setwd("C:\\Users\\Deng\\Dropbox\\MI R codes")
#' Multiple Imputation through Indirect Use of Regularized Regression(IURR) for univariate missing data pattern
#' 
#' @param data  A data frame or a matrix containing the incomplete data. Missing values are coded as NA.
#' @param family Specify the type of the missing data. Can be "gaussian", or "binary", or "poisson".
#' @param method The method to be used in regularized regression. Can be "lasso", "elastic net","adaptive lasso" or "blasso".
#' @param m Number of multiple imputations. The default is m=5.
#' @return n.missing Number of missing values in the data.
#' @return col.missing Which column contains the missing values.
#' @return impute The imputations(Column number equals to m). Each column corresponds to one imputation.
#' @export
#' @examples
#' data(datagaussian)
#' data(databinary)
#' data(datapoisson)
#' MIiurr(data=datagaussian,method="lasso",family="gaussian")
#' MIiurr(data=datagaussian,method="blasso",family="gaussian")
#' MIiurr(data=databinary,method="elastic net",family="binary")
#' MIiurr(data=databinary,method="adaptive lasso",family="binary")
#' MIiurr(data=datapoisson,method="elastic net",family="poisson")
#' MIiurr(data=datapoisson,method="adaptive lasso",family="poisson")



MIiurr=function(data,family=c("gaussian"),method=c("lasso"),m=5)
{
  set.seed(123)
  realdata=as.matrix(data)
  n=nrow(realdata)
  
  colflag=as.numeric(which(colSums(is.na(realdata))!=0))
  x1=realdata[,colflag]
  nmis=sum(is.na(x1))
  newobs<-realdata[,-colflag]
  inpute=matrix(NA,n,m)
  
  res<-realdata[,colflag]
  obs<-realdata[,-colflag]
  res<-as.matrix(res)
  obs<-as.matrix(obs)
  
  res1=res[!is.na(x1),]
  obs1=obs[!is.na(x1),]
  
  obs2=obs[is.na(x1),]
  #######################################################
  #Gaussian##############################################
  #######################################################
  
  if(family=="gaussian"){
    ##lasso
    if(method=="lasso"){
      nonz=c(1)
      while(length(nonz)==1 | length(nonz)>=40){
      model<- glmnet(obs1,res1,standardize=FALSE)
      cvres<-cv.glmnet(obs1,res1,standardize=FALSE)
      fits<-coef(model, s=cvres$lambda.min)
      fits<-as.numeric(fits)
      #print(is.logical(fits!=0))
      nonz<-which(fits!=0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-lm(res1~obs1[,nonzero]),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-newcoff[!is.na(newcoff)]
      
      newcov<-vcov(newreg)
      nonzero<-as.data.frame(nonzero)
      
    }
    
    ##elastic
    if(method=="elastic net"){
      
      model<- glmnet(obs1,res1,standardize=FALSE,alpha=0.5)
      cvres<-cv.glmnet(obs1,res1,standardize=FALSE,alpha=0.5)
      fits<-coef(model, s=cvres$lambda.min)
      nonz<-which(fits!=0)
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-lm(res1~obs1[,nonzero]),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-na.omit(newcoff)
      newcov<-vcov(newreg)
      nonzero<-as.data.frame(nonzero)
      
    }
    
    ##adaptive lasso
    if(method=="adaptive lasso"){
      
      model<-adalasso(obs1,res1)
      
      bbb<-model$coefficients.adalasso
      nonzero<-as.data.frame(which(bbb!=0))
      
      try(newreg<-lm(res1~obs1[,nonzero[[1]]]),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-na.omit(newcoff)
      newcov<-vcov(newreg)
    }
    
    ###blasso
    
    if(method=="blasso"){
      model<-blasso.vs(res1,obs1,sig2=1,tau=1,phi=0.5,sig2prior=c(0.1,0.1),thin=2, tauprior=c(0.01,0.01),iters=10000,beta=rep(0.1,ncol(realdata)))
      bbb<-model$beta[sample(5000:nrow(model$beta),m,replace=F),]             
    }
    if(method=="blasso"){
      mu4<-0
      for(k in 1:m){
        nonzero<-which(bbb[k,]!=0)          
        if(length(nonzero)>1){
          inter<-mean(res1)-apply(obs1[,nonzero],2,mean)%*%(bbb[k,which(bbb[k,]!=0)]/apply(obs1[,nonzero],2,sd))
          predd<-rep(inter,length(res1))+obs1[,nonzero]%*%(bbb[k,which(bbb[k,]!=0)]/apply(obs1[,nonzero],2,sd))
          sse<-sqrt(sum((res1-predd)^2)/n)
          
          pred<-rep(inter,nrow(obs2))+obs2[,nonzero]%*%(bbb[k,which(bbb[k,]!=0)]/apply(obs2[,nonzero],2,sd))
        }else{
          inter<-mean(res1)-mean(obs1[,nonzero])*(bbb[k,which(bbb[k,]!=0)]/sd(obs1[,nonzero]))
          predd<-rep(inter,length(res1))+obs1[,nonzero]*(bbb[k,which(bbb[k,]!=0)]/sd(obs1[,nonzero]))
          sse<-sqrt(sum((res1-predd)^2)/n)
          pred<-rep(inter,nrow(obs2))+obs2[,nonzero]*(bbb[k,which(bbb[k,]!=0)]/sd(obs2[,nonzero]))          
        }
        inpute[,k]<-x1
        ii=0
        for(i in 1:n){
          if(is.na(inpute[i,k])){
            ii=ii+1
            mu4<-pred[ii]
            inpute[i,k]<-rnorm(1,mu4,sd=sse)
          }
        }
      }       
    }else{
      ##inpute missing values
      mu4<-0
      for(k in 1:m){
        newmodel<-mvrnorm(n = 1, newcoff, newcov, tol = 1e-6, empirical = FALSE)
        if(length(newmodel)>2){
          predd<-newmodel[1]+obs1[,nonzero[1:length(newmodel)-1,]]%*%newmodel[2:length(newmodel)]
          sse<-sqrt(sum((res1-predd)^2)/length(res1))
          
          pred<-newmodel[1]+obs2[,nonzero[1:length(newmodel)-1,]]%*%newmodel[2:length(newmodel)]
        }else{
          predd<-newmodel[1]+obs1[,nonzero[[1]]]*newmodel[2]
          sse<-sqrt(sum((res1-predd)^2)/length(res1))
          pred<-newmodel[1]+obs2[,nonzero[[1]]]*newmodel[2]
        }
        
        inpute[,k]<-x1
        ii=0
        for(i in 1:n){
          if(is.na(inpute[i,k])){
            ii=ii+1
            mu4<-pred[ii]
            inpute[i,k]<-rnorm(1,mu4,sd=sse)
          }
        }
      }
      
    }
  }   
  
  #######################################################
  #Binary##############################################
  #######################################################
  if(family=="binary"){
    ##lasso
    if(method=="lasso"){
      nonz=1
      while(length(nonz)==1){
        model<- glmnet(obs1,res1,family="binomial",standardize=FALSE)
        cvres<-cv.glmnet(obs1,res1,family="binomial",standardize=FALSE)
        fits<-coef(model, s=cvres$lambda.min)
        nonz<-which(fits!=0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-glm(res1~obs1[,nonzero],family="binomial"),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-na.omit(newcoff)
      newcov<-vcov(newreg)
      nonzero<-as.data.frame(nonzero)
      
    }
    
    ##elastic
    if(method=="elastic net"){
      nonz=1
      while(length(nonz)==1){
        model<- glmnet(obs1,res1,family="binomial",standardize=FALSE,alpha=0.5)
        cvres<-cv.glmnet(obs1,res1,family="binomial",standardize=FALSE,alpha=0.5)
        fits<-coef(model, s=cvres$lambda.min)
        fits_vec <- as.numeric(fits)         # This flattens fits to a numeric vector
        nonz <- which(fits_vec != 0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-glm(res1~obs1[,nonzero],family="binomial"),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-na.omit(newcoff)
      newcov<-vcov(newreg)
      nonzero<-as.data.frame(nonzero)
      
    }
    
    ##adaptive lasso
    if(method=="adaptive lasso"){
      nonz=1
      while(length(nonz)==1){
        model <- glmnet(obs1, res1, family="binomial", standardize=FALSE, alpha=0)
        betals <- coef(model)
        betals_vec <- as.numeric(betals)
        penalty.factor <- 1 / abs(betals_vec[-1])
        penalty.factor[!is.finite(penalty.factor)] <- 1
        
        if(length(penalty.factor) != ncol(obs1)) {
          warning("Forcing penalty.factor to match ncol(obs1)")
          penalty.factor <- penalty.factor[seq_len(ncol(obs1))]
        }
        cvres <- cv.glmnet(obs1, res1, family="binomial", standardize=FALSE, alpha=1, penalty.factor=penalty.factor)
        fits<-coef(model, s=cvres$lambda.min)
        fits_vec <- as.numeric(fits)         # This flattens fits to a numeric vector
        nonz <- which(fits_vec != 0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg <- glm(res1 ~ obs1[, nonzero], family="binomial"), silent=TRUE)
      idx <- which(!is.na(newreg$coefficients))
      newcoff <- newreg$coefficients[idx]
      newcov <- vcov(newreg)[idx, idx, drop=FALSE]
      nonzero <- as.data.frame(nonzero)
    }        
    
    ##inpute missing values
    mu4<-0
    for(k in 1:m){
      newmodel<-mvrnorm(n = 1, newcoff, newcov, tol = 1e-6, empirical = FALSE)
      if(length(newmodel)>2){
        # predd<-newmodel[1]+obs1[,nonzero[1:length(newmodel)-1,]]%*%newmodel[2:length(newmodel)]
        pred<-newmodel[1]+obs2[,nonzero[1:length(newmodel)-1,]]%*%newmodel[2:length(newmodel)]
      }else{
        #predd<-newmodel[1]+obs1[,nonzero[[1]]]*newmodel[2]
        pred<-newmodel[1]+obs2[,nonzero[[1]]]*newmodel[2]
      }
      
      inpute[,k]<-x1
      ii=0
      for(i in 1:n){
        if(is.na(inpute[i,k])){
          ii=ii+1
          if(pred[ii]>100){mu4=1
          }else{
            mu4<-exp(pred[ii])/(1+exp(pred[ii]))
          }
          inpute[i,k]<-rbinom(1,1,mu4)
        }
      }
    }
    
  }
  
  
  
  #######################################################
  #Poisson##############################################
  #######################################################
  if(family=="poisson"){
    ##lasso
    if(method=="lasso"){
      nonz=1
      while(length(nonz)==1){
        model<- glmnet(obs1,res1,family="poisson",standardize=FALSE)
        cvres<-cv.glmnet(obs1,res1,family="poisson",standardize=FALSE)
        fits<-coef(model, s=cvres$lambda.min)
        nonz<-which(fits!=0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-glm(res1~obs1[,nonzero],family="poisson"),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-na.omit(newcoff)
      newcov<-vcov(newreg)
      nonzero<-as.data.frame(nonzero)
      
    }
    
    ##elastic
    if(method=="elastic net"){
      nonz=1
      while(length(nonz)==1){
        model<- glmnet(obs1,res1,family="poisson",standardize=FALSE,alpha=0.5)
        cvres<-cv.glmnet(obs1,res1,family="poisson",standardize=FALSE,alpha=0.5)
        fits<-coef(model, s=cvres$lambda.min)
        fits_vec <- as.numeric(fits)         # This flattens fits to a numeric vector
        nonz <- which(fits_vec != 0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-glm(res1~obs1[,nonzero],family="poisson"),silent=TRUE)
      newcoff<-newreg$coefficients
      newcoff<-na.omit(newcoff)
      newcov<-vcov(newreg)
      nonzero<-as.data.frame(nonzero)
      
    }
    
    ##adaptive lasso
    if(method=="adaptive lasso"){
      nonz=1
      while(length(nonz)==1){
        model<- glmnet(obs1,res1,family="poisson",standardize=FALSE,alpha=0)
        betals=coef(model)
        betals_vec <- as.numeric(betals)
        penalty.factor <- 1 / abs(betals_vec[-1])
        penalty.factor[!is.finite(penalty.factor)] <- 1
        
        if(length(penalty.factor) != ncol(obs1)) {
          warning("Forcing penalty.factor to match ncol(obs1)")
          penalty.factor <- penalty.factor[seq_len(ncol(obs1))]
        }
        cvres<-cv.glmnet(obs1,res1,family="poisson",standardize=FALSE,alpha=1,penalty.factor=penalty.factor,)
        fits<-coef(model, s=cvres$lambda.min)
        fits_vec <- as.numeric(fits)         # This flattens fits to a numeric vector
        nonz <- which(fits_vec != 0)
      }
      nonzero<-nonz[-1]-1
      #if(length(nonzero)>length(res1)-10){next}
      try(newreg<-glm(res1~obs1[,nonzero],family="poisson"),silent=TRUE)
      # Subset to non-NA coefficients
      idx <- which(!is.na(newreg$coefficients))
      newcoff <- newreg$coefficients[idx]
      newcov <- vcov(newreg)[idx, idx, drop = FALSE]
      nonzero<-as.data.frame(nonzero)
    }        
    
    ##inpute missing values
    mu4<-0
    for(k in 1:m){
      newmodel<-mvrnorm(n = 1, newcoff, newcov, tol = 1e-6, empirical = FALSE)
      if(length(newmodel)>2){
        #predd<-newmodel[1]+obs1[,nonzero[1:length(newmodel)-1,]]%*%newmodel[2:length(newmodel)]
        pred<-newmodel[1]+obs2[,nonzero[1:length(newmodel)-1,]]%*%newmodel[2:length(newmodel)]
      }else{
        #predd<-newmodel[1]+obs1[,nonzero[[1]]]*newmodel[2]
        pred<-newmodel[1]+obs2[,nonzero[[1]]]*newmodel[2]
      }
      
      inpute[,k]<-x1
      ii=0
      for(i in 1:n){
        if(is.na(inpute[i,k])){
          ii=ii+1
          mu4<-exp(pred[ii])
          inpute[i,k]<-rpois(1,mu4)
        }
      }
    }
    
  }    
  
  
  
  
  relist=list(n.missing=nmis,col.missing=colflag,impute=inpute)
  
  
  return(relist)
  
  
}#end of function