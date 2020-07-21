##############
#Data Reading#
##############

#setwd("workdirectory")
library(caret)
library(data.table)
library(glmnet)
library(ModelMetrics)
library(tidyr)
library(tree)

test.r=fread("test_V2.csv")
train.r=fread("train_V2.csv")

#basic info
head(test.r)
dim(test.r) #1934174 28
head(train.r)
dim(train.r) #4446966 29

################
#Data Cleansing#
################

#there exist some players with 0 in ride distance and walk distance,there should be a bug
cheater=train.r[train.r$rideDistance==0 & train.r$walkDistance==0 & 
                  train.r$swimDistance==0 & train.r$kills!=0,] 
dim(cheater) #1535 29
train.nocheater=train.r[-as.numeric(row.names(cheater)),]
summary(train.nocheater[,-c(1,2,3)])
str(train.nocheater)

#remove some no used variables
train=train.nocheater[,-c(1,2,3)]

#change -1 in train$rankPoints as 0
train[train$rankPoints==-1,16]=0

#check if na inside and remove them
train[!complete.cases(train),]
train=train[complete.cases(train),]

#save train data
save(train,file="train.rda")

################
#Data Exploring#
################

#load train data
load("train.rda")

#Distribution of WinPlacePerc
hist(train$winPlacePerc,main="Distribution of WinPlacePerc",xlab="WinPlacePerc",
     col="antiquewhite4",breaks=40)
box()

#Pearson Correlation
cor.m=cor(train[,-13])
cor.m[lower.tri(cor.m)]=NA
mcor.m=melt(cor.m,na.rm=TRUE)
ggheatmap=ggplot(data=mcor.m, aes(Var2, Var1, fill=value))+
  geom_tile(color="white") + 
  scale_fill_gradient2(low="antiquewhite4", high="forestgreen", mid="white", 
                       midpoint=0, limit=c(-1,1), space="Lab", name="Pearson Correlation") + 
  theme_minimal() + theme(axis.text.x=element_text(angle=45, vjust=1, size=8, hjust=1)) + 
  coord_fixed() + ggtitle("PUBG Correlation Plot")
ggheatmap + geom_text(aes(Var2, Var1, label=round(value, 2)), color="black", size=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_blank(), 
        panel.border=element_blank(), panel.background=element_blank(), axis.ticks=element_blank(), 
        legend.justification=c(1, 0), legend.position=c(0.6, 0.7), legend.direction="horizontal")+
  guides(fill=guide_colorbar(barwidth=7, barheight=1,title.position="top", title.hjust=0.5))

################
#Data Splitting#
################

#table for matchType
table(train$matchType)

#split train data according to matchtype
crashfpp=train[grep("^crashfpp$",train$matchType),-13] #6285 25
crashtpp=train[grep("^crashtpp$",train$matchType),-13] #371  25
duo=train[grep("^duo$",train$matchType),-13]
duofpp=train[grep("^duo-fpp$",train$matchType),-13]
flarefpp=train[grep("^flarefpp$",train$matchType),-13]
flaretpp=train[grep("^flaretpp$",train$matchType),-13]
normalduo=train[grep("^normal-duo$",train$matchType),-13]
normalduofpp=train[grep("^normal-duo-fpp$",train$matchType),-13]
normalsolo=train[grep("^normal-solo$",train$matchType),-13]
normalsolofpp=train[grep("^normal-solo-fpp$",train$matchType),-13]
normalsquad=train[grep("^normal-squad$",train$matchType),-13]
normalsquadfpp=train[grep("^normal-squad-fpp$",train$matchType),-13]
solo=train[grep("^solo$",train$matchType),-13]
solofpp=train[grep("^solo-fpp$",train$matchType),-13]
squad=train[grep("^squad$",train$matchType),-13]
squadfpp=train[grep("^squad-fpp$",train$matchType),-13]

################
#Model Choosing#
################

#function for linear
li=function(x.train,x.test,y.train,y.test){
  model.l=vector("list",ncol(y.train))
  x.train=as.data.frame(x.train)
  x.test=as.data.frame(x.test)
  y.train=as.data.frame(y.train)
  pred.m=as.data.frame(y.test)
  coef.m=matrix(NA,ncol(x.train)+1,ncol(y.train))
  for(ii in 1:ncol(y.test)){
    model.l[[ii]]=lm(y.train[,ii]~.,data=x.train)
    coef.m[,ii]=as.numeric(coef(model.l[[ii]]))
    pred.m[,ii]=predict(model.l[[ii]],newdata=x.test)
    pred.m[pred.m[,ii]<0,ii]=0
    pred.m[pred.m[,ii]>1,ii]=1
  }
  return(list(pred.m=pred.m,model.l=model.l,coef.m=coef.m))
}

#function for glmnet with choosing alpha 
gl=function(x.train,x.test,y.train,y.test){
  alphas=c(0:10)/10
  pred.m=as.data.frame(y.test)
  x.train=as.matrix(x.train)
  x.test=as.matrix(x.test)
  y.train=as.data.frame(y.train)
  model.l=vector("list",ncol(y.train))
  ealpha=matrix(NA,nrow=1,ncol=ncol(y.train))
  coef.m=matrix(NA,ncol(x.train)+1,ncol(y.train))
  for(ii in 1:ncol(y.train)){
    tmp0=vector("list",length(alphas))
    foldid = createFolds(y.train[,ii], k=5, list=F)
    for(aa in 1:length(alphas)){
      tmp0[[aa]]=cv.glmnet(x=x.train,y=y.train[,ii], foldid=foldid, alpha=alphas[aa],family="gaussian")
    } 
    ord=which.min(sapply(tmp0, function(x.train) x.train$cvm[which(x.train$lambda==x.train$lambda.1se)]))
    ealpha[1,ii]=alphas[ord]
    model.l[[ii]]=tmp0[[ord]]
    coef.m[,ii]=as.numeric(coef(model.l[[ii]],s="lambda.1se"))
    pred.m[,ii]=predict(model.l[[ii]],newx=x.test,s="lambda.1se")
    pred.m[pred.m[,ii]<0,ii]=0
    pred.m[pred.m[,ii]>1,ii]=1
    }
    return(list(pred.m=pred.m, model.l=model.l, coef.m=coef.m, alpha=ealpha))
}

#function for regression tree
tr=function(x.train,x.test,y.train,y.test){
  model.l=vector("list",ncol(y.train))
  pred.m=as.data.frame(y.test)
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
    pred.m[,ii]=predict(model.l[[ii]],x.test)
    pred.m[pred.m[,ii]<0,ii]=0
    pred.m[pred.m[,ii]>1,ii]=1
  }
  return(list(pred.m=pred.m,model.l=model.l))
}

#function for Mean Absolute Error
maer=function(method,data,nfold){
  yy=data[,25,drop=FALSE]
  xx=data[,-25]
  pred=as.data.frame(yy)
  coef.m=matrix(NA,ncol(xx)+1,nfold)
  model.l=vector("list",nfold)
  idx.cv=createFolds(as.matrix(yy[,1]),k=nfold,list=F)
  for(ll in 1:nfold){
    y.train = yy[idx.cv!=ll,] 
    y.test = yy[idx.cv==ll, ]   
    x.train = xx[idx.cv!=ll, ]
    x.test = xx[idx.cv==ll, ]
    if(method=="linear"){
      bundle=li(x.train,x.test,y.train,y.test)
      coef.m[,ll]=bundle$coef.m
    } else if(method=="glmnet"){
      bundle=gl(x.train,x.test,y.train,y.test)
      coef.m[,ll]=bundle$coef.m
    } else if(method=="regressiontree"){
      bundle=tr(x.train,x.test,y.train,y.test)
    }
    pred[idx.cv==ll,]=bundle$pred.m
    model.l[[ll]]=bundle$model.l[[1]]
  }
  mae.n=mae(as.matrix(yy),as.matrix(pred))
  return(list(mae.n=mae.n,coef.m=coef.m,model.l=model.l))
}

#function for finding Mean Absolute Error for different matchType 
nmaer=function(data){
	time1=Sys.time()
	if(nrow(data)<=500){
	  nrun=100
	  nf=5
	} else if (nrow(data)<=10000){
	  nrun=100
	  nf=10
	} else if(nrow(data)<=100000){
	  nrun=20
	  nf=10
	} else if(nrow(data)>100000){
	  nrun=20
	  nf=5
	}
	cat(nrow(data),nrun,nf,"\n")
	mae.m=matrix(NA,nrow=nrun,ncol=3)
	set.seed(0)
	zz=sample(1:1000,nrun)
	rr=0
	for(ss in zz){
  		set.seed(ss)
  		rr=rr+1
  		mae.m[rr,1]=maer("linear",data,nf)$mae.n
  		cat(rr,1,"finish","\n")
  		mae.m[rr,2]=maer("glmnet",data,nf)$mae.n
  		cat(rr,2,"finish","\n")
  		mae.m[rr,3]=maer("regressiontree",data,nf)$mae.n
  		cat(rr,3,"finish","\n")
	}
	mae.m=as.data.frame(mae.m)
	names(mae.m)=c("linear","glmnet","regressiontree")
	row.names(mae.m)=zz
	time2=Sys.time()
	print(time2-time1)
	return(mae.m)
}

#get each match type's Mean Absolute Error and save data
crashfpp.mae=nmaer(crashfpp) #Time difference of 33.21087 mins
crashtpp.mae=nmaer(crashtpp) #Time difference of 4.851793 mins
duo.mae=nmaer(duo) #Time difference of 3.525404 hours
duofpp.mae=nmaer(duofpp) #Time difference of 10.30849 hours
flarefpp.mae=nmaer(flarefpp) #Time difference of 16.85413 mins
flaretpp.mae=nmaer(flaretpp) #Time difference of 21.79857 mins
normalduo.mae=nmaer(normalduo) #Time difference of 7.050116 mins
normalduofpp.mae=nmaer(normalduofpp) #Time difference of 35.43435 mins
normalsolo.mae=nmaer(normalsolo) #Time difference of 6.339052 mins
normalsolofpp.mae=nmaer(normalsolofpp) #Time difference of 14.98273 mins
normalsquad.mae=nmaer(normalsquad) #Time difference of 17.76672 mins
normalsquadfpp.mae=nmaer(normalsquadfpp) #Time difference of 17.4742 mins
solo.mae=nmaer(solo) #Time difference of 1.501289 hours
solofpp.mae=nmaer(solofpp) #Time difference of 4.352886 hours
squad.mae=nmaer(squad) #Time difference of 5.537296 hours
squadfpp.mae=nmaer(squadfpp) #Time difference of 17.11544 hours

save(crashfpp.mae,crashtpp.mae,duo.mae,duofpp.mae, flarefpp.mae,flaretpp.mae,normalduo.mae,normalduofpp.mae,normalsolo.mae,normalsolofpp.mae,normalsquad.mae,normalsquadfpp.mae,solo.mae,solofpp.mae,squad.mae,squadfpp.mae,file="mae.rda")

#load data
load("mae.rda")

#function for choosing method and model 
cho=function(data.mae,data,name){
  nc=order(colMeans(data.mae))[1]
  minseed=as.numeric(rownames(data.mae[order(data.mae[,nc])[1],]))
  set.seed(minseed)
  if(nrow(data)<=500){
    nfold=5
    nrun=100
  }else if (nrow(data)<=10000){
    nfold=10
    nrun=100
  } else if(nrow(data)<=100000){
    nfold=10
    nrun=20
  } else if(nrow(data)>100000){
    nfold=5
    nrun=20
  }
  boxplot(data.mae,main=name)
  legend("topleft",legend=c(paste("nrun=",nrun),paste("nfold=",nfold)))
  if(names(data.mae)[nc] %in% c("linear","glmnet")){
    coef.m=rowMeans(maer(names(data.mae)[nc],data,nfold)$coef.m)
    coef.m[is.na(coef.m)]=0
    return(list(method=names(data.mae)[nc],coef.m=coef.m))
  }else if(names(data.mae)[nc] %in% c("regressiontree")){
    model.l= maer(names(data.mae)[nc],data,nfold)$model.l
    pred.m=matrix(NA,nrow(data),nfold)
    mae.n=matrix(NA,nrow=nfold,ncol=1)
    for(jj in 1:nfold){
      yy=data[,25,drop=FALSE]
      xx=data[,-25]
      pred.m[,jj]=predict(model.l[[jj]],xx)
      mae.n[jj,]=mae(as.matrix(yy),pred.m[,jj])
    }
    model.b=model.l[[order(mae.n[,1])[1]]]
    return(list(method=names(data.mae)[nc],model.b=model.b))
  }
}

#choose method and model for each match type
crash1=cho(crashfpp.mae,crashfpp,"crashfpp.mae")
crash3=cho(crashtpp.mae,crashtpp,"crashtpp.mae")
duo3=cho(duo.mae,duo,"duotpp.mae")
duo1=cho(duofpp.mae,duofpp,"duofpp.mae")
flare1=cho(flarefpp.mae,flarefpp,"flarefpp.mae")
flare3=cho(flaretpp.mae,flaretpp,"flaretpp.mae")
nduo3=cho(normalduo.mae,normalduo,"normalduotpp.mae")
nduo1=cho(normalduofpp.mae,normalduofpp,"normalduofpp.mae")
nsolo3=cho(normalsolo.mae,normalsolo,"normalsolotpp.mae")
nsolo1=cho(normalsolofpp.mae,normalsolofpp,"normalsolofpp.mae")
nsquad3=cho(normalsquad.mae,normalsquad,"nnormalsquadtpp.mae")
nsquad1=cho(normalsquadfpp.mae,normalsquadfpp,"nnormalsquadfpp.mae")
solo3=cho(solo.mae,solo,"solotpp.mae")
solo1=cho(solofpp.mae,solofpp,"solofpp.mae")
squad3=cho(squad.mae,squad,"squadtpp.mae")
squad1=cho(squadfpp.mae,squadfpp,"squadfpp.mae")

save(crash1,crash3,duo3,duo1,flare1,flare3,nduo3,nduo1,nsolo3,nsolo1,nsquad3,nsquad1,solo3,solo1,squad3,squad1, file="models.rda")

################
#Data Analyzing#
################

#function for color choosing
color=function(data){
  col=NA
  for(cc in 1:length(data)){
    if(data[cc]>=0){
      col=c(col,"forestgreen")
    } else if(data[cc]<0){
      col=c(col,"antiquewhite4")
    }
  }
  return(col[-1])
}
#barplot for each match type's coefficients
barplot(crash1$coef.m[-1],main="Coefficients of crash-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals",
                    "killPlace","killPoints","kills","killStreaks","longestKill","matchDuration",
                    "maxPlace","numGroups","rankPoints","revives","rideDistance","roadKills",
                    "swimDistance","teamKills","vehicleDestroys","walkDistance","weaponsAcquired",
                    "winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(crash1$coef.m[-1]))
barplot(crash3$coef.m[-1],main="Coefficients of crash-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(crash3$coef.m[-1]))
barplot(duo3$coef.m[-1],main="Coefficients of duo-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(duo3$coef.m[-1]))
barplot(duo1$coef.m[-1],main="Coefficients of duo-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(duo1$coef.m[-1]))
barplot(flare1$coef.m[-1],main="Coefficients of flare-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(flare1$coef.m[-1]))
barplot(flare3$coef.m[-1],main="Coefficients of flare-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(flare3$coef.m[-1]))
barplot(nduo3$coef.m[-1],main="Coefficients of normal-duo-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(nduo3$coef.m[-1]))
barplot(nduo1$coef.m[-1],main="Coefficients of normal-duo-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(nduo1$coef.m[-1]))
barplot(nsolo3$coef.m[-1],main="Coefficients of normal-solo-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(nsolo3$coef.m[-1]))
plot(nsolo1$model.b,col="antiquewhite4")
text(nsolo1$model.b,pretty=0,cex=0.8,col="forestgreen")
title("Coefficients of normal-solo-fpp Plot")
barplot(nsquad3$coef.m[-1],main="Coefficients of normal-squad-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(nsquad3$coef.m[-1]))
barplot(nsquad1$coef.m[-1],main="Coefficients of normal-squad-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(nsquad1$coef.m[-1]))
barplot(solo3$coef.m[-1],main="Coefficients of normal-solo-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(solo3$coef.m[-1]))
barplot(solo1$coef.m[-1],main="Coefficients of normal-solo-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(solo1$coef.m[-1]))
barplot(squad3$coef.m[-1],main="Coefficients of squad-tpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(squad3$coef.m[-1]))
barplot(squad1$coef.m[-1],main="Coefficients of squad-fpp Plot",
        names.arg=c("assists","boosts","damageDealt","DBNOs","headshotKills","heals","killPlace",
                    "killPoints","kills","killStreaks","longestKill","matchDuration","maxPlace",
                    "numGroups","rankPoints","revives","rideDistance","roadKills","swimDistance",
                    "teamKills","vehicleDestroys","walkDistance","weaponsAcquired","winPoints"),
        horiz=FALSE,las=2,cex.names=0.5,yaxp=c(-0.15,0.05,4),col=color(squad1$coef.m[-1]))



#################
#Data Predicting#
#################

#load models
load("models.rda")

#function for predicting each match type's value
predf=function(model,rr){
  if(model$method %in% c("linear","glmnet")){
    pred=as.matrix(cbind(1,test[rr,-13])) %*% as.matrix(model$coef.m)
  } else if (model$method %in% c("regressiontree")){
    pred=predict(model$model.b,test[rr,-13])
  }
  return(pred)
}

#remove some no used variables
test=test.r[,-c(1,2,3)]
summary(test)

#change -1 in test$rankPoints as 0
test[test$rankPoints==-1,16]=0

#use different methods to predict test data
pred=matrix(NA,nrow=nrow(test),ncol=1)
time1=Sys.time()
for(pp in 1:nrow(test)){
  if(test[pp,13]=="crashfpp") pred[pp,1]=predf(crash1,pp)
  if(test[pp,13]=="crashtpp") pred[pp,1]=predf(crash3,pp)
  if(test[pp,13]=="duo") pred[pp,1]=predf(duo3,pp)
  if(test[pp,13]=="duo-fpp") pred[pp,1]=predf(duo1,pp)
  if(test[pp,13]=="flarefpp") pred[pp,1]=predf(flare1,pp)
  if(test[pp,13]=="flaretpp") pred[pp,1]=predf(flare3,pp)
  if(test[pp,13]=="normal-duo") pred[pp,1]=predf(nduo3,pp)
  if(test[pp,13]=="normal-duo-fpp") pred[pp,1]=predf(nduo1,pp)
  if(test[pp,13]=="normal-solo") pred[pp,1]=predf(nsolo3,pp)
  if(test[pp,13]=="normal-solo-fpp") pred[pp,1]=predf(nsolo1,pp)
  if(test[pp,13]=="normal-squad") pred[pp,1]=predf(nsquad3,pp)
  if(test[pp,13]=="normal-squad-fpp") pred[pp,1]=predf(nsquad1,pp)
  if(test[pp,13]=="solo") pred[pp,1]=predf(solo3,pp)
  if(test[pp,13]=="solo-fpp") pred[pp,1]=predf(solo1,pp)
  if(test[pp,13]=="squad") pred[pp,1]=predf(squad3,pp)
  if(test[pp,13]=="squad-fpp") pred[pp,1]=predf(squad1,pp)
}
time2=Sys.time()
print(time2-time1) #Time difference of 4.937633 hours

#check pred
pred[is.na(pred)]
pred[pred[,1]<0,1]=0
pred[pred[,1]>1,1]=1

#store results
write.csv(cbind(test.r$Id,pred),file="sample_submission_V2.csv")






