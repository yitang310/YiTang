#setwd and call library
setwd("/Users/ty/Desktop/stat454/hw/hw2/A2_2_Tang_Yi_(and_Yang_Chenxi)")
source("HelperFunctions.R")
library(caret)
library(class)
library(glmnet)
library(MASS)
library(tidyr)
library(tree)
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

#xx for dataframe all mutation, yy for matrix for ic50, yy is log10 for yyn.
xx=data.frame(data.3tc[,2:229])
yy=data.frame(data.3tc$Y,data.abc$Y,data.azt$Y,data.d4t$Y,data.ddi$Y)
yyn=10^yy

#dataframe for cutoff point
cutoff=data.frame(3,2,3,1.5,1.5)
names(cutoff)=c("3tc","abc","azt","d4t","ddi")

#yynew is the dataframe for categorical variables
cc3tc=factor(yyn$data.3tc.Y>cutoff[,1],labels=c("susceptible","resistance"))
ccabc=factor(yyn$data.abc.Y>cutoff[,2],labels=c("susceptible","resistance"))
ccazt=factor(yyn$data.azt.Y>cutoff[,3],labels=c("susceptible","resistance"))
ccd4t=factor(yyn$data.d4t.Y>cutoff[,4],labels=c("susceptible","resistance"))
ccddi=factor(yyn$data.ddi.Y>cutoff[,5],labels=c("susceptible","resistance"))
yynew=data.frame(cc3tc,ccabc,ccazt,ccd4t,ccddi)

#function for the proportion cutoff
pcut=function(pred,y.test,minn,maxx){
  min.n=ceiling(minn*100)
  max.x=floor(maxx*100)
  tmp0=matrix(NA,nrow=(max.x-min.n+1),ncol=1)
  kk=0
  for(jj in min.n:max.x){
    kk=kk+1
    tmp1=factor(pred>jj/100,labels=c("susceptible","resistance"))
    tmp0[kk,1]=mean(tmp1!=y.test)
  }
  ord=(order(tmp0[,1])[1]+min.n-1)/100
  return(ord)
}

#function for logistic regression
lo=function(x.train,x.test,y.train,y.test){
  pred.m=matrix(NA,nrow=nrow(y.test),ncol=ncol(y.test))
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  pp=matrix(NA,nrow=1,ncol=ncol(y.train))
  model.l=vector("list",ncol(y.train))
  for(ii in 1:ncol(y.test)){
    pred.train=matrix(NA,nrow=nrow(y.train),ncol=1)
    model.l[[ii]]=glm(y.train[,ii]~.,family="binomial",data=x.train)
    pred.train[,1]=predict(model.l[[ii]],x.train,type="response")
    pp[1,ii]=pcut(pred.train,y.train[,ii],min(pred.train[,1]),max(pred.train[,1]))
    pred.m[,ii]=predict(model.l[[ii]],x.test,type="response")
  }
  return(list(pred.m=pred.m,model.l=model.l,ppcut=pp))
}

#function for logistic glmnet
logl=function(x.train,x.test,y.train,y.test,alpha){
  pred.m=matrix(NA,nrow=nrow(y.test),ncol=ncol(y.test))
  x.train=as.matrix(x.train)
  x.test=as.matrix(x.test)
  pp=matrix(NA,nrow=1,ncol=ncol(y.train))
  model.l=vector("list",ncol(y.train))
  for(ii in 1:ncol(y.test)){
    pred.train=matrix(NA,nrow=nrow(y.train),ncol=1)
    model.l[[ii]]=cv.glmnet(x.train, y.train[,ii], alpha=alpha,family="binomial") 
    pred.train[,1]=predict(model.l[[ii]],x.train,type="response")
    pp[1,ii]=pcut(pred.train,y.train[,ii],min(pred.train[,1]),max(pred.train[,1]))
    pred.m[,ii]=predict(model.l[[ii]],x.test,type="response")
  }
  return(list(pred.m=pred.m,model.l=model.l,ppcut=pp))
} 

#function for logistic glmnet with choosing alpha 
logl.better=function(x.train,x.test,y.train,y.test){
  alphas=c(0:10)/10
  pp=matrix(NA,nrow=1,ncol=ncol(y.train))
  pred.m=matrix(NA,nrow=nrow(y.test),ncol=ncol(y.test))
  x.train=as.matrix(x.train)
  x.test=as.matrix(x.test)
  model.l=vector("list",ncol(y.train))
  ealpha=matrix(NA,nrow=1,ncol=ncol(y.train))
  for(ii in 1:ncol(y.train)){
    tmp0=vector("list",length(alphas))
    pred.train=matrix(NA,nrow=nrow(y.train),ncol=1)
    foldid = createFolds(y.train[,ii], k=5, list=F)
    for(aa in 1:length(alphas)){
      tmp0[[aa]]=cv.glmnet(x=x.train,y=y.train[,ii], foldid=foldid, alpha=alphas[aa],family="binomial")
    } 
    ord=which.min(sapply(tmp0, function(x.train) x.train$cvm[which(x.train$lambda==x.train$lambda.1se)]))
    ealpha[1,ii]=alphas[ord]
    model.l[[ii]]=tmp0[[ord]]
    pred.train[,1]=predict(model.l[[ii]],x.train,type="response")
    pp[1,ii]=pcut(pred.train,y.train[,ii],min(pred.train[,1]),max(pred.train[,1]))
    pred.m[,ii]=predict(model.l[[ii]],x.test,type="response")
    
  }
  return(list(pred.m=pred.m,model.l=model.l,each.alpha=ealpha,ppcut=pp))
}

#function for Linear Discriminant Analysis
ldaa=function(x.train,x.test,y.train,y.test){
  pred.class=as.data.frame(y.test)
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  model.l=vector("list",ncol(y.train))
  for(ii in 1:ncol(y.test)){
    model.l[[ii]]=lda(y.train[,ii]~.,data=x.train) 
    pred=predict(model.l[[ii]],x.test)
    pred.class[,ii]=pred$class
  }
  return(list(pred.m=pred.class,model.l=model.l))
}

#function for choosing K-Nearest Neighbors's k
kch=function(x.train,y.train){
  kk=1+2*c(0:9)
  tmp0=matrix(NA,nrow=length(kk),ncol=1)
  qq=0
  for(jj in kk){
    qq=qq+1
    tmp1=knn(x.train,x.train,y.train,k=jj)
    tmp0[qq,1]=mean(tmp1!=y.train)
  }
  ord=(order(tmp0[,1])[1]-1)*2+1
  return(ord)
}

#function for K-Nearest Neighbors
knnn=function(x.train,x.test,y.train,y.test){
  model.l=vector("list",ncol(y.train)) #not exist
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  y.train=as.data.frame(y.train)
  pred.m=as.data.frame(y.test)
  for(ii in 1:ncol(y.test)){
    k=kch(x.train,y.train[,ii])
    pred.m[,ii]=knn(x.train,x.test,y.train[,ii],k=k)
  }
  return(list(pred.m=pred.m,model.l=model.l))
}

#function for classification tree
cltr=function(x.train,x.test,y.train,y.test){
  model.l=vector("list",ncol(y.train))
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  y.train=as.data.frame(y.train)
  pred.m=as.data.frame(y.test)
  assign("x.train", x.train, .GlobalEnv)
  assign("y.train",y.train, .GlobalEnv)
  for(ii in 1:ncol(y.test)){
    assign("ii", ii, .GlobalEnv)
    tmp0=tree(y.train[,ii]~.,data = x.train)
    cv.tmp0=cv.tree(tmp0,FUN=prune.misclass) #prune the tree
    nn=cv.tmp0$size[order(cv.tmp0$dev)[1]]
    if(nn==1) nn=2 #if best equal 1, it will cause "singlenode" problem in predict
    model.l[[ii]]=prune.misclass(tmp0,best=nn)  
    pred.m[,ii]=predict(model.l[[ii]],x.test,type="class")
  }
  return(list(pred.m=pred.m,model.l=model.l))
}

#function for linear
#use log value to predict first and exponentiate back will make classification rate better
li=function(x.train,x.test,y.train,y.test){
  model.l=vector("list",ncol(y.train))
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  y.train=as.data.frame(y.train)
  pred.m=as.data.frame(y.test)
  for(ii in 1:ncol(y.test)){
    model.l[[ii]]=lm(y.train[,ii]~.,data=x.train)
    pred=10^(predict(model.l[[ii]],newdata=x.test))
    pred.m[,ii]=rep("susceptible",nrow(y.test))
    pred.m[pred>cutoff[,ii],ii]="resistance"
  }
  return(list(pred.m=pred.m,model.l=model.l))
}

#function for glmnet
#use log value to predict first and exponentiate back will make classification rate better
gl=function(x.train,x.test,y.train,y.test,alpha){
  pred.m=as.data.frame(y.test)
  x.train=as.matrix(x.train)
  x.test=as.matrix(x.test)
  model.l=vector("list",ncol(y.train))
  for(ii in 1:ncol(y.test)){
    model.l[[ii]]=cv.glmnet(x.train, y.train[,ii], alpha=alpha,family="gaussian") 
    pred=10^(predict(model.l[[ii]],newx=x.test,s="lambda.1se"))
    pred.m[,ii]=rep("susceptible",nrow(y.test))
    pred.m[pred>cutoff[,ii],ii]="resistance"
  }
  return(list(pred.m=pred.m,model.l=model.l))
} 

#function for glmnet with choosing alpha 
#use log value to predict first and exponentiate back will make classification rate higher
gl.better=function(x.train,x.test,y.train,y.test){
  alphas=c(0:10)/10
  pred.m=as.data.frame(y.test)
  x.train=as.matrix(x.train)
  x.test=as.matrix(x.test)
  model.l=vector("list",ncol(y.train))
  ealpha=matrix(NA,nrow=1,ncol=ncol(y.train))
  for(ii in 1:ncol(y.train)){
    tmp0=vector("list",length(alphas))
    foldid = createFolds(y.train[,ii], k=5, list=F)
    for(aa in 1:length(alphas)){
      tmp0[[aa]]=cv.glmnet(x=x.train,y=y.train[,ii], foldid=foldid, alpha=alphas[aa],family="gaussian")
    } 
    ord=which.min(sapply(tmp0, function(x.train) x.train$cvm[which(x.train$lambda==x.train$lambda.1se)]))
    ealpha[1,ii]=alphas[ord]
    model.l[[ii]]=tmp0[[ord]]
    pred=10^(predict(model.l[[ii]],newx=x.test,s="lambda.1se"))
    pred.m[,ii]=rep("susceptible",nrow(y.test))
    pred.m[pred>cutoff[,ii],ii]="resistance"
    
  }
  return(list(pred.m=pred.m,model.l=model.l,alpha=ealpha))
}

#function for regression tree
tr=function(x.train,x.test,y.train,y.test){
  model.l=vector("list",ncol(y.train))
  pred.m=as.matrix(y.test)
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  y.train=as.data.frame(y.train)
  assign("x.train", x.train, .GlobalEnv)
  assign("y.train",y.train, .GlobalEnv)
  for(ii in 1:ncol(y.test)){
    assign("ii", ii, .GlobalEnv)
    tmp0=tree(y.train[,ii]~.,data=x.train)
    cv.tmp0=cv.tree(tmp0)   #prune the tree
    nn=cv.tmp0$size[order(cv.tmp0$dev)[1]]
    if(nn==1) nn=2 #if best equal 1, it will cause "singlenode" problem in predict
    model.l[[ii]]=prune.tree(tmp0,best=nn) 
    pred=10^(predict(model.l[[ii]],x.test))
    pred.m[,ii]=rep("susceptible",nrow(y.test))
    pred.m[pred>cutoff[,ii],ii]="resistance"
  }
  return(list(pred.m=pred.m,model.l=model.l))
}

#function for all methods
nmethod=function(method,x.train,x.test,y.train,y.test,alpha=NULL){
  ppcut=NA
  if(method=="logistic") {
    fit=lo(x.train,x.test,y.train,y.test)
    pred.m=fit$pred.m
    ppcut=fit$ppcut
  }
  if(method == "logisticglmnet"){
    if(is.null(alpha)){
      fit=logl.better(x.train,x.test,y.train,y.test)
      pred.m=fit$pred.m
      ppcut=fit$ppcut
    } else{
      fit=logl(x.train,x.test,y.train,y.test,alpha)
      pred.m=fit$pred.m
      ppcut=fit$ppcut
    }
  }
  if(method=="lda"){
    pred.m=ldaa(x.train,x.test,y.train,y.test)$pred.m
  }
  if(method=="knn"){
    pred.m=knnn(x.train,x.test,y.train,y.test)$pred.m
  }
  if(method=="classificationtree"){
    pred.m=cltr(x.train,x.test,y.train,y.test)$pred.m
  }
  if(method=="linear") {
    pred.m=li(x.train,x.test,y.train,y.test)$pred.m
  }
  if(method=="glmnet"){
    if(is.null(alpha)){
      pred.m=gl.better(x.train,x.test,y.train,y.test)$pred.m
    }else{
      pred.m=gl(x.train,x.test,y.train,y.test,alpha)$pred.m
    }
  }
  if(method=="regressiontree") {
    pred.m=tr(x.train,x.test,y.train,y.test)$pred.m
  }
  return(list(pred.m=pred.m,ppcut=ppcut))
}

#function for misclassification rate
misclassification=function(method,nfold,alpha=NULL){
  if(method %in% c("logistic","logisticglmnet","lda","knn","classificationtree")){
    x.m=xx
    y.m=yynew
    if(method %in% c("logistic","logisticglmnet")){
      pred.m=matrix(NA,nrow=nrow(y.m),ncol=ncol(y.m))
    } else {
      pred.m=as.data.frame(y.m)
    }
    
  } else if(method %in% c("linear","glmnet","regressiontree")){
    x.m=xx
    y.m=yy
    pred.m=as.data.frame(y.m)
  }
  idx.cv=createFolds(y.m[,1],k=nfold,list=F)
  ppcut=matrix(NA, nrow=nfold,ncol=ncol(y.m))
  for(ll in 1:nfold){
    y.train = y.m[idx.cv!=ll,]   
    y.test = y.m[idx.cv==ll, ]   
    x.train = x.m[idx.cv!=ll, ]
    x.test = x.m[idx.cv==ll, ]
    tmp=nmethod(method,x.train,x.test,y.train,y.test,alpha)
    ppcut[ll,]=tmp$ppcut
    pred.m[idx.cv==ll,]=tmp$pred.m
  }
  if(method %in% c("logistic","logisticglmnet")){
    pred=pred.m
    pred.m=matrix(NA,nrow=nrow(y.m),ncol=ncol(y.m))
    for(jj in 1:ncol(y.m)){
      pp=mean(ppcut[,jj])
      pred.m[,jj]=rep("susceptible",nrow(y.m))
      pred.m[pred[,jj]>pp,jj]="resistance"
    }
  }
  misrate=mean(pred.m!=yynew)
  return(misrate)
}

#since cv.glmnet runs very slow, I set a example seed to find the best alpha in logisticglmnet
##compare logistic ridge, logistic lasso and logistic elastic net
set.seed(204)
misclassification("logisticglmnet",10,alpha=0) #0.1044944
misclassification("logisticglmnet",10,alpha=1) #0.09775281
misclassification("logisticglmnet",10,alpha=0.1) #0.09951846
misclassification("logisticglmnet",10,alpha=0.2) #0.09566613
misclassification("logisticglmnet",10,alpha=0.3) #0.09357945
misclassification("logisticglmnet",10,alpha=0.4) #0.09919743
misclassification("logisticglmnet",10,alpha=0.5) #0.09486356
misclassification("logisticglmnet",10,alpha=0.6) #0.09357945
misclassification("logisticglmnet",10,alpha=0.7) #0.09839486
misclassification("logisticglmnet",10,alpha=0.8) #0.09855538
misclassification("logisticglmnet",10,alpha=0.9) #0.09550562
set.seed(620)
misclassification("logisticglmnet",10,alpha=0.3) #0.09438202
misclassification("logisticglmnet",10,alpha=0.6) #0.09582665

#run 100 times to find best method which is giving lowest misclassification rate
zz=c(1:100)
misclassificationrate=matrix(NA,nrow=100,ncol=8)
for(ss in zz){
  set.seed(ss)
  misclassificationrate[ss,1]=misclassification("logistic",10)
  cat(ss,1,"finish","\n")
  misclassificationrate[ss,2]=misclassification("logisticglmnet",10,alpha=0.3)
  #misclassificationrate[ss,2]=misclassification("logisticglmnet",10)
  cat(ss,2,"finish","\n")
  misclassificationrate[ss,3]=misclassification("lda",10)
  cat(ss,3,"finish","\n")
  misclassificationrate[ss,4]=misclassification("knn",10)
  cat(ss,4,"finish","\n")
  misclassificationrate[ss,5]=misclassification("classificationtree",10)
  cat(ss,5,"finish","\n")
  misclassificationrate[ss,6]=misclassification("linear",10)
  cat(ss,6,"finish","\n")
  #misclassificationrate[ss,7]=misclassification("glmnet",10,alpha=0.1)
  misclassificationrate[ss,7]=misclassification("glmnet",10)
  cat(ss,7,"finish","\n")
  misclassificationrate[ss,8]=misclassification("regressiontree",10)
  cat(ss,8,"finish","\n")
}
misclassificationrate=as.data.frame(misclassificationrate)
names(misclassificationrate)=c("logistic","logisticglmnet 0.3 ","lda","knn","classificationtree","linear","glmnet","regressiontree")
save(misclassificationrate, file="misclassificationrate.rda")

#boxplot
load("misclassificationrate.rda")
boxplot(misclassificationrate[1:8])


