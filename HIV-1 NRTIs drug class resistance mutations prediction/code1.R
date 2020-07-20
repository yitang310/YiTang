#setwd and call library
setwd("/Users/ty/Desktop/stat454/hw/hw1/A1_2_Tang_Yi_(and_Yang_Chenxi)")
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion")
source("HelperFunctions.R")
library(caret)
library(glmnet)
library(lava)
library(tidyr)

#read data
mut.TSM = mut.name(Mut="TSM.NRTI")
mut.expert = mut.name(Mut="Exp.NRTI")
mut.comp = mut.name(Mut="Comp.NRTI")
data_3tc = load.data(dataset="NRTI",drug = "3TC",min.muts = 10,muts.in = mut.comp)
data_abc = load.data(dataset="NRTI",drug = "ABC",min.muts = 10,muts.in = mut.comp)
data_azt = load.data(dataset="NRTI",drug = "AZT",min.muts = 10,muts.in = mut.comp)
data_d4t = load.data(dataset="NRTI",drug = "D4T",min.muts = 10,muts.in = mut.comp)
data_ddi = load.data(dataset="NRTI",drug = "DDI",min.muts = 10,muts.in = mut.comp)

#Getting common row name and column name 
##Function for getting common row.names
newrow=function(a,b,c,d,e){   
  ord=numeric(0)
  for (i in 1:nrow(a)){
    if((row.names(a)[i] %in% row.names(b) == TRUE) & 
       (row.names(a)[i] %in% row.names(c) == TRUE) & 
       (row.names(a)[i] %in% row.names(d) == TRUE) & 
       (row.names(a)[i] %in% row.names(e) == TRUE)){
      ord=append(ord,row.names(a)[i])
    } else{
      cat("Row", row.names(a)[i], "excluded from the model because it do not appear in other drugs' data.\n")
    }
  }
  return(ord)
}

##Getting the common row.names
new.row=newrow(data_3tc,data_abc,data_azt,data_d4t,data_ddi)  

##Function for getting a new drug's data that each drug's data has same row names
data_newdrugy=function(drug,rownumber){  
  return(drug[rownumber,])
}

##Getting 5 new drugs' data which all have same row.names.
data_3tcy=data_newdrugy(data_3tc,new.row)
data_abcy=data_newdrugy(data_abc,new.row)
data_azty=data_newdrugy(data_azt,new.row)
data_d4ty=data_newdrugy(data_d4t,new.row)
data_ddiy=data_newdrugy(data_ddi,new.row)

##Function for getting common column names
newcol=function(aa,bb,cc,dd,ee){   
  ord=numeric(0)
  for (ii in 1:ncol(aa)){
    if((names(aa)[ii] %in% names(bb) == TRUE) & 
       (names(aa)[ii] %in% names(cc) == TRUE) & 
       (names(aa)[ii] %in% names(dd) == TRUE) & 
       (names(aa)[ii] %in% names(ee) == TRUE)){
      ord=append(ord,names(aa)[ii])
    }else{
      cat(names(aa)[ii], "excluded from the model because it do not appear in other drugs' data.")
    }
  }
  return(ord)
}

##Getting the common column names
new.col=newcol(data_3tcy,data_abcy,data_azty,data_d4ty,data_ddiy) 

##Function for getting a new drug's data that each drug's data has same column names
data_newdrug=function(drug,colnumber){  
  return(drug[,colnumber])
}

##Getting 5 new drugs' data which all have same column names.
##(same column names and same row names)
data.3tc=data_newdrug(data_3tcy,new.col)
data.abc=data_newdrug(data_abcy,new.col)
data.azt=data_newdrug(data_azty,new.col)
data.d4t=data_newdrug(data_d4ty,new.col)
data.ddi=data_newdrug(data_ddiy,new.col)

#xx for dataframe all mutation, yy for matrix for ic50
xx=data.frame(data.3tc[,2:229])
yy=data.frame(data.3tc$Y,data.abc$Y,data.azt$Y,data.d4t$Y,data.ddi$Y)

#function for linear
li=function(x.m,y.m){
  cy=ncol(y.m)
  x.m=as.data.frame(x.m)
  yhat.m<-matrix(NA,nrow=nrow(y.m),ncol=cy)
  coef.m<-matrix(NA,nrow=ncol(x.m)+1,ncol=cy)
  for(ii in 1:cy){
    ll=lm(y.m[,ii]~.,data=x.m)
    coef.m[,ii]=as.numeric(coef(ll))
    yhat.m[,ii]=predict(ll,x.m)
  }
  return(list(yhat.m=yhat.m, coef.m=coef.m))
}

#function for glmnet
gl=function(x.m,y.m,alpha){
  cy=ncol(y.m)
  yhat.m<-matrix(NA,nrow=nrow(y.m),ncol=cy)
  coef.m<-matrix(NA,ncol(x.m)+1,cy)
  for(ii in 1:cy){
    gg=cv.glmnet(as.matrix(x.m),y.m[,ii],alpha=alpha)
    coef.m[,ii]=as.numeric(coef(gg,s="lambda.1se"))
    yhat.m[,ii]=predict(gg,newx=as.matrix(x.m),s="lambda.1se")
  }
  return(list(yhat.m=yhat.m, coef.m=coef.m))
}

#function for no cv 
no.cv=function(method,x.m,y.m,alpha=NULL){
  if (method == "linear"){
    coef.m=li(x.m,y.m)$coef.m
    yhat.m=li(x.m,y.m)$yhat.m
  } else if (method %in% c("ridge","lasso","elasticnet")){
    coef.m=gl(x.m,y.m,alpha)$coef.m
    yhat.m=gl(x.m,y.m,alpha)$yhat.m
  }
  return(list(yhat.m=yhat.m, coef.m=coef.m))
}

#function for cv linear
cv.li=function(x.m,y.m,nfold){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  yhat.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  coef.a=array(NA,dim=c(ncol(x.m)+1,cy,nfold))
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]
    y.test = y.m[idx.cv==jj, ]
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    coef.m = li(x.train,y.train)$coef.m
    yhat.m[idx.cv==jj,] = as.matrix(cbind(1,x.test)) %*% coef.m
    coef.a[,,jj] = coef.m
  }
  return(list(yhat.m=yhat.m, coef.m=apply(coef.a, c(1,2),mean) ))
}

#function for cv glmnet
cv.gl=function(x.m,y.m,alpha,nfold){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  yhat.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  coef.a=array(NA,dim=c(ncol(x.m)+1,cy,nfold))
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]
    y.test = y.m[idx.cv==jj, ]
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    coef.m = gl(x.train,y.train,alpha=alpha)$coef.m
    yhat.m[idx.cv==jj,] = as.matrix(cbind(1,x.test)) %*% coef.m
    coef.a[,,jj] = coef.m
  }
  return(list(yhat.m=yhat.m, coef.m=apply(coef.a, c(1,2),mean) ))
}

# function for have cv
yes.cv=function(method,x.m,y.m,alpha=NULL,nfold=NULL){
  if (method == "linear"){
    coef.m=cv.li(x.m,y.m,nfold)$coef.m
    yhat.m=cv.li(x.m,y.m,nfold)$yhat.m
  } else if (method %in% c("ridge","lasso","elasticnet")){
    coef.m=cv.gl(x.m,y.m,alpha,nfold)$coef.m
    yhat.m=cv.gl(x.m,y.m,alpha,nfold)$yhat.m
  }
  return(list(yhat.m=yhat.m, coef.m=coef.m))
}

##########################
#function for no stacking#
##########################
no.stack=function(method,x.m,y.m,alpha=NULL,nfold=NULL){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  yhat.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  error.m=yhat.m
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]   
    y.test = y.m[idx.cv==jj, ]   
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    coef.m=no.cv(method,x.train,y.train,alpha=alpha)$coef.m
    yhat.m[idx.cv==jj,]=as.matrix(cbind(1,x.test)) %*% coef.m
  }
  error.m=yhat.m-as.matrix(y.m)
  mse=mean(error.m^2)
  variance=tr(var(error.m))
  bias=mean(error.m)
  is.na(alpha)=NA
  result=result=data.frame(method,alpha,mse,variance,bias)
  names(result)=c("model","alpha","average MSE","total variance","average bias")
  row.names(result)="no stacking"
  return(result)
}

################################
#function for standard stacking#
################################
st.stack=function(method,x.m,y.m,alpha=NULL,nfold=NULL){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  ytilde.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  error.m=ytilde.m
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]   
    y.test = y.m[idx.cv==jj, ]   
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    nocv=no.cv(method,x.train,y.train,alpha=alpha)
    coef.beta=nocv$coef.m
    coef.alpha=no.cv(method,nocv$yhat.m,nocv$yhat.m,alpha=alpha)$coef.m
    ytilde.m[idx.cv==jj,]=as.matrix(cbind(1,as.matrix(cbind(1,x.test)) %*% coef.beta)) %*% coef.alpha
  }
  error.m=ytilde.m-as.matrix(y.m)
  mse=mean(error.m^2)
  variance=tr(var(error.m))
  bias=mean(error.m)
  is.na(alpha)=NA
  result=result=data.frame(method,alpha,mse,variance,bias)
  names(result)=c("model","alpha","average MSE","total variance","average bias")
  row.names(result)="standard stacking"
  return(result)
}

#########################################################
#function for standard stacking with linear in 2nd stage#
#########################################################
stl.stack=function(method,x.m,y.m,alpha=NULL,nfold=NULL){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  ytilde.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  error.m=ytilde.m
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]   
    y.test = y.m[idx.cv==jj, ]   
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    nocv=no.cv(method,x.train,y.train,alpha=alpha)
    coef.beta=nocv$coef.m
    coef.alpha=no.cv("linear",nocv$yhat.m,nocv$yhat.m,alpha=alpha)$coef.m
    ytilde.m[idx.cv==jj,]=as.matrix(cbind(1,as.matrix(cbind(1,x.test)) %*% coef.beta)) %*% coef.alpha
  }
  error.m=ytilde.m-as.matrix(y.m)
  mse=mean(error.m^2)
  variance=tr(var(error.m))
  bias=mean(error.m)
  is.na(alpha)=NA
  result=result=data.frame(method,"linear",alpha,mse,variance,bias)
  names(result)=c("model1","model2","alpha","average MSE","total variance","average bias")
  row.names(result)="standard stacking"
  return(result)
}

###################################
#function for cv improved stacking#
###################################
cv.stack=function(method,x.m,y.m,alpha=NULL,nfold=NULL){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  ytilde.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  error.m=ytilde.m
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]   
    y.test = y.m[idx.cv==jj, ]   
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    yescv=yes.cv(method,x.train,y.train,alpha=alpha,nfold=nfold)
    coef.beta=yescv$coef.m
    coef.alpha=no.cv(method,yescv$yhat.m,yescv$yhat.m,alpha=alpha)$coef.m
    ytilde.m[idx.cv==jj,]=as.matrix(cbind(1,as.matrix(cbind(1,x.test)) %*% coef.beta)) %*% coef.alpha
  }
  error.m=ytilde.m-as.matrix(y.m)
  mse=mean(error.m^2)
  variance=tr(var(error.m))
  bias=mean(error.m)
  is.na(alpha)=NA
  result=result=data.frame(method,alpha,mse,variance,bias)
  names(result)=c("model","alpha","average MSE","total variance","average bias")
  row.names(result)="cv improved stacking"
  return(result)
}

############################################################
#function for cv improved stacking with linear in 2nd stage#
############################################################
cvl.stack=function(method,x.m,y.m,alpha=NULL,nfold=NULL){
  cy=ncol(y.m)
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  ytilde.m=matrix(NA,nrow=nrow(y.m),ncol=cy)
  error.m=ytilde.m
  for(jj in 1:nfold){
    y.train = y.m[idx.cv!=jj,]   
    y.test = y.m[idx.cv==jj, ]   
    x.train = x.m[idx.cv!=jj, ]
    x.test = x.m[idx.cv==jj, ]
    yescv=yes.cv(method,x.train,y.train,alpha=alpha,nfold=nfold)
    coef.beta=yescv$coef.m
    coef.alpha=no.cv("linear",yescv$yhat.m,yescv$yhat.m,alpha=alpha)$coef.m
    ytilde.m[idx.cv==jj,]=as.matrix(cbind(1,as.matrix(cbind(1,x.test)) %*% coef.beta)) %*% coef.alpha
  }
  error.m=ytilde.m-as.matrix(y.m)
  mse=mean(error.m^2)
  variance=tr(var(error.m))
  bias=mean(error.m)
  is.na(alpha)=NA
  result=result=data.frame(method,"linear",alpha,mse,variance,bias)
  names(result)=c("model1","model2","alpha","average MSE","total variance","average bias")
  row.names(result)="cv improved stacking"
  return(result)
}

#########
#Results#
#########
#no stacking
no.stack("linear",xx,yy,nfold=10)
no.stack("ridge",xx,yy,alpha=0,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.1,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.2,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.3,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.4,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.5,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.6,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.7,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.8,nfold=10)
no.stack("elasticnet",xx,yy,alpha=0.9,nfold=10)
no.stack("lasso",xx,yy,alpha=1,nfold=10)

#standard stacking
st.stack("linear",xx,yy,nfold=10)
st.stack("ridge",xx,yy,alpha=0,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.1,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.2,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.3,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.4,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.5,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.6,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.7,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.8,nfold=10)
st.stack("elasticnet",xx,yy,alpha=0.9,nfold=10)
st.stack("lasso",xx,yy,alpha=1,nfold=10)

#cv improved stacking
cv.stack("linear",xx,yy,nfold=10)
cv.stack("ridge",xx,yy,alpha=0,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.1,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.2,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.3,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.4,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.5,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.6,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.7,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.8,nfold=10)
cv.stack("elasticnet",xx,yy,alpha=0.9,nfold=10)
cv.stack("lasso",xx,yy,alpha=1,nfold=10)

#standard stacking with linear in 2nd stage
stl.stack("linear",xx,yy,nfold=10)
stl.stack("ridge",xx,yy,alpha=0,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.1,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.2,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.3,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.4,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.5,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.6,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.7,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.8,nfold=10)
stl.stack("elasticnet",xx,yy,alpha=0.9,nfold=10)
stl.stack("lasso",xx,yy,alpha=1,nfold=10)

#cv improved stacking with linear in 2nd stage
cvl.stack("linear",xx,yy,nfold=10)
cvl.stack("ridge",xx,yy,alpha=0,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.1,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.2,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.3,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.4,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.5,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.6,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.7,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.8,nfold=10)
cvl.stack("elasticnet",xx,yy,alpha=0.9,nfold=10)
cvl.stack("lasso",xx,yy,alpha=1,nfold=10)