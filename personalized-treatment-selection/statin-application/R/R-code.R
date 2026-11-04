#################################################################
library(MASS)
library(nnet)
library(glmnet)
library(twang)
library(gbm)
library(tibble)
library(dplyr)
library(kedd)
library(fastDummies)


##Function to obtain adaptive LASSO model and the selected variables
xnames<-function(df){
  x=model.matrix(y~.,df[,-ncol(df)])[,-1]
  y=df$y
  ridge <- glmnet(x = x, y = y,alpha = 0)
  ## Perform ridge regression with 10-fold CV
  ridge_cv <- cv.glmnet(x = x, y = y,type.measure = "mse",nfold = 10,alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.min))[-1]
  
  ##Adaptive LASSO
  alasso <- glmnet(x = x,y = y,alpha=1,penalty.factor = 1/abs(best_ridge_coef))
  alasso_cv <- cv.glmnet(x = x, y = y,type.measure = "mse",nfold = 10,alpha = 1,
                         penalty.factor = 1 / abs(best_ridge_coef),keep = TRUE)
  best_alasso_coef <- coef(alasso_cv, s = alasso_cv$lambda.min)
  lam<-alasso_cv$lambda.min
  #las = function(x0) { predict(alasso_cv, s=lam, newx = x0 ) }
  keep.x<-names(which(abs(best_alasso_coef[,1][-1])>0))
  return(list(keep.x=keep.x,alasso_cv=alasso_cv))
}


#Function to tune the parameters in the GBM and return the final model and selected variables based on relative influence
var.gbm<-function(df){
  
  # create grid search
  hyper_grid1 <- expand.grid(
    learning_rate = c(0.3, 0.1, 0.05, 0.01, 0.005),
    RMSE = NA,
    trees = NA,
    time = NA
  )
  
  # execute grid search
  for(i in seq_len(nrow(hyper_grid1))) {
    
    # fit gbm
    train_time <- system.time({
      mdl <- gbm(
        formula = y ~ .,
        data = subset(df, select=-c(z)),
        distribution = "gaussian",
        n.trees = 5000, 
        shrinkage = hyper_grid1$learning_rate[i], 
        interaction.depth = 3, 
        n.minobsinnode = 10,
        cv.folds = 5
      )
    })
    
    # add SSE, trees, and training time to results
    hyper_grid1$RMSE[i]  <- sqrt(min(mdl$cv.error))
    hyper_grid1$trees[i] <- which.min(mdl$cv.error)
    hyper_grid1$Time[i]  <- train_time[["elapsed"]]
    
  }
  
  ind<-which.min(hyper_grid1$RMSE)
  
  mdl <- gbm(
    formula = y ~ .,
    data = subset(df, select=-c(z)),
    distribution = "gaussian",
    n.trees = hyper_grid1$trees[ind], 
    shrinkage = hyper_grid1$learning_rate[ind], 
    interaction.depth = 3, 
    n.minobsinnode = 10,
    cv.folds = 5
  )
  
  rel.inf<-relative.influence(mdl, hyper_grid1$trees[ind], scale. = TRUE, sort. = TRUE)/sum(relative.influence(mdl, hyper_grid1$trees[ind], scale. = TRUE, sort. = TRUE))
  
  keep.x<-names(rel.inf[rel.inf>0.05])
  return(list(trees=hyper_grid1$trees[ind],
              shrinkage = hyper_grid1$learning_rate[ind],
              keep.x=keep.x,mdl=mdl))
}

#Local Linear Estimator and NW estimator-
alpha.hat<-function(dff,si,di,zi,h){
  dat<-subset(dff,z==zi & d==di,select=c(y,S,d,z))
  if (nrow(dat)<5){
    return(list(nwe=-100,lle=-100))
  }else{
    y<-dat$y
    S<-dat$S
    d<-dat$d
    z<-dat$z
    #h<-h.amise(S,deriv.order = 1)$h
    Z<-cbind(1,S-si)
    nwe<-(((1/h)*dnorm((S-si)/h))%*%y)/(sum((1/h)*dnorm((S-si)/h)))
    zWz<-t(Z)%*%diag((1/h)*dnorm((S-si)/h))%*%Z
    
    if (det(zWz)<1e-5) {
      lle<--100
    }else{
      lle<-c(1,0)%*%solve(t(Z)%*%diag((1/h)*dnorm((S-si)/h))%*%Z)%*%(t(Z)%*%diag((1/h)*dnorm((S-si)/h))%*%y)
    }
    return(list(nwe=nwe,lle=lle))
  }}


#Local Linear Estimator and NW estimator-weighted
alpha.hat.w<-function(dff,si,di,zi,h){
  dat<-subset(dff,z==zi & d==di,select=c(y,S,d,w,z))
  if (nrow(dat)<5){
    return(list(nwe=-100,lle=-100))
  }else{
    y<-dat$y
    S<-dat$S
    d<-dat$d
    w<-dat$w
    z<-dat$z
    Z<-cbind(1,S-si)
    
    nwe<-(((1/h)*w*dnorm((S-si)/h))%*%y)/(sum((1/h)*w*dnorm((S-si)/h)))
    zWz<-t(Z)%*%diag(w*(1/h)*dnorm((S-si)/h))%*%Z
    
    if (det(zWz)<1e-5) {
      lle<--100
    }else{
      lle<-c(1,0)%*%solve(t(Z)%*%diag(w*(1/h)*dnorm((S-si)/h))%*%Z)%*%(t(Z)%*%diag(w*(1/h)*dnorm((S-si)/h))%*%y)
    }
    return(list(nwe=nwe,lle=lle))
  }}


#Method that selects
match.gbm<-function(df,xn){
  df1=subset(df,z==1)
  df2=subset(df,z==2)
  df3=subset(df,z==3)
  #df3=subset(df,z==3)
  # keep.x1<-sel.gbm(df1,5000,0.01)
  # keep.x2<-sel.gbm(df2,5000,0.01)
  # keep.x3<-sel.gbm(df3,5000,0.01)
  # keep.x4<-sel.gbm(df4,5000,0.01)
  
  gbm.sel1<-var.gbm(df1)
  keep.x1<-gbm.sel1$keep.x
  best.iter10<-gbm.sel1$trees
  shrink1<-gbm.sel1$shrinkage
  mdl1<-gbm.sel1$mdl
  #best.iter10<- gbm.perf(mdl1, method = "cv",plot.it = FALSE)
  
  gbm.sel2<-var.gbm(df2)
  keep.x2<-gbm.sel2$keep.x
  best.iter20<-gbm.sel2$trees
  shrink2<-gbm.sel2$shrinkage
  mdl2<-gbm.sel2$mdl
  #best.iter20<- gbm.perf(mdl2, method = "cv",plot.it = FALSE)
  
  gbm.sel3<-var.gbm(df3)
  keep.x3<-gbm.sel3$keep.x
  best.iter30<-gbm.sel3$trees
  shrink3<-gbm.sel3$shrinkage
  mdl3<-gbm.sel3$mdl
  #best.iter30<- gbm.perf(mdl3, method = "cv",plot.it = FALSE)
  
  #keep.x3<-xnames(df3)
  keep.x<-unique(c(keep.x1,keep.x2,keep.x3))
  #ntrees<-5000
  #shrink<-0.01
  
  x.all<-names(subset(df, select=-c(y,z)))
  
  drop.x<-x.all[!x.all%in%keep.x]
  dfw<-df[ , !(names(df) %in% drop.x)]
  
  dff<-df
  dfw$z<-relevel(dfw$z,ref=3)
  m1<<-multinom(z~.,data = dfw[,-1],trace=FALSE)
  pr<-predict(m1,type = "probs")
  #dff$p1<-predict(m1,type = "probs")[,1]
  #dff$p2<-predict(m1,type = "probs")[,2]
  #dff$p3<-1-dff$p1-dff$p2
  pr.z<-numeric(nrow(dfw))
  for (i in 1:nrow(dfw)) {
    pr.z[i]<-pr[i,dff$z[i]]
    
  }
  
  psFormula <- paste(keep.x, collapse="+")
  psFormula <- formula(paste("z~",psFormula, sep=""))
  mnps.gbm <- mnps(psFormula, data=dfw[,-1],
                   estimand="ATE",n.trees=5000,stop.method="es.max",verbose = FALSE)
  
  summary(mnps.gbm)
  dff$w <- get.weights(mnps.gbm, estimand = "ATE",
                       stop.method = "es.max")
  
  #Scores
  Sa<-predict(mdl1,newdata = dff,n.trees =best.iter10)-apply(cbind(predict(mdl2,newdata = dff,n.trees =best.iter20),predict(mdl3,newdata = dff,n.trees =best.iter30)),1,max)
  Sb<-predict(mdl2,newdata = dff,n.trees =best.iter20)-apply(cbind(predict(mdl1,newdata = dff,n.trees =best.iter10),predict(mdl3,newdata = dff,n.trees =best.iter30)),1,max)
  Sc<-predict(mdl3,newdata = dff,n.trees =best.iter30)-apply(cbind(predict(mdl1,newdata = dff,n.trees =best.iter10),predict(mdl2,newdata = dff,n.trees =best.iter20)),1,max)
  
  dff$S<-apply(cbind(Sa,Sb,Sc),1,max)
  dff$d<-factor(apply(cbind(Sa,Sb,Sc),1,which.max),levels=c(1,2,3))
  
  ind1<-which(table(subset(dff,z==1)$d)>4)
  ind2<-which(table(subset(dff,z==2)$d)>4)
  ind3<-which(table(subset(dff,z==3)$d)>4)
  
  dff1<-subset(subset(dff,z==1),d%in%c(ind1))
  dff2<-subset(subset(dff,z==2),d%in%c(ind2))
  dff3<-subset(subset(dff,z==3),d%in%c(ind3))
  
  h1<-h.amise(dff1$S,deriv.order = 1)$h
  h2<-h.amise(dff2$S,deriv.order = 1)$h
  h3<-h.amise(dff3$S,deriv.order = 1)$h
  
  #Score model
  #my<-lm(y~S*d,data = dff)
  my1<-lm(y~S*d,weights = w,data = dff1)
  my2<-lm(y~S*d,weights = w,data = dff2)
  my3<-lm(y~S*d,weights = w,data = dff3)
  
  S0a<-predict(mdl1,newdata = xn,n.trees =best.iter10)-apply(cbind(predict(mdl2,newdata = xn,n.trees =best.iter20),predict(mdl3,newdata = xn,n.trees =best.iter30)),1,max)
  S0b<-predict(mdl2,newdata = xn,n.trees =best.iter20)-apply(cbind(predict(mdl1,newdata = xn,n.trees =best.iter10),predict(mdl3,newdata = xn,n.trees =best.iter30)),1,max)
  S0c<-predict(mdl3,newdata = xn,n.trees =best.iter30)-apply(cbind(predict(mdl1,newdata = xn,n.trees =best.iter10),predict(mdl2,newdata = xn,n.trees =best.iter20)),1,max)
  S0<-apply(cbind(S0a,S0b,S0c),1,max)
  d0<-factor(apply(cbind(S0a,S0b,S0c),1,which.max),levels = c(1,2,3))
  
  check<-data.frame(matrix(nrow=m,ncol = 3))
  tt<-data.frame(matrix(nrow=m,ncol = 2))
  nw.lam<-data.frame(matrix(nrow=m,ncol = 2))
  ll.lam<-data.frame(matrix(nrow=m,ncol = 2))
  nwest<-data.frame(matrix(nrow=m,ncol = 4))
  for(k in 1:m){
    
    # nonpar1<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=1,h1)
    # nonpar2<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=2,h2)
    # nonpar3<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=3,h3)
    nonpar1<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=1,h1)
    nonpar2<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=2,h2)
    nonpar3<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=3,h3)
    
    ynw1<-nonpar1$nwe
    ynw2<-nonpar2$nwe
    ynw3<-nonpar3$nwe
    
    yll1<-nonpar1$lle
    yll2<-nonpar2$lle
    yll3<-nonpar3$lle
    
    ynw<-sort(c(ynw1,ynw2,ynw3))
    nw.lam[k,]<-c(ynw[3]-ynw[2],ynw[3]-ynw[1])
    
    yll<-sort(c(yll1,yll2,yll3))
    ll.lam[k,]<-c(yll[3]-yll[2],yll[3]-yll[1])
    
    nwest[k,]<-c(ynw1,ynw2,ynw3,which.max(c(ynw1,ynw2,ynw3)))
    check[k,]<-c(yn[k],which.max(c(ynw1,ynw2,ynw3)),which.max(c(yll1,yll2,yll3)))
    tt[k,1]<-ifelse(check[k,1]==check[k,2],1,0)
    tt[k,2]<-ifelse(check[k,1]==check[k,3],1,0)
  }
  C.nw<-max(nw.lam)
  C.ll<-max(ll.lam)
  
  nw.Lam12<-apply(nw.lam,2,mean)/C.nw
  ll.Lam12<-apply(ll.lam,2,mean)/C.ll
  nw.assign<-table(yn,check[,2])
  ll.assign<-table(yn,check[,3])
  accprob<-apply(tt,2,mean)
  return(list(nw.assign=nw.assign,ll.assign=ll.assign,accprob=accprob,nw.Lam12=nw.Lam12,ll.Lam12=ll.Lam12))
}


match.lasso<-function(df,xn){
  x.all<-names(subset(df, select=-c(y,z)))
  df1=subset(df,z==1)
  df2=subset(df,z==2)
  df3=subset(df,z==3)
  
  obj1<-xnames(df1)
  las1<-obj1$alasso_cv
  
  obj2<-xnames(df2)
  las2<-obj2$alasso_cv
  
  obj3<-xnames(df3)
  las3<-obj3$alasso_cv
  
  keep.x1<-obj1$keep.x
  keep.x2<-obj2$keep.x
  keep.x3<-obj3$keep.x
  
  keep.x<-unique(c(keep.x1,keep.x2,keep.x3))
  
  drop.x<-x.all[!x.all%in%keep.x]
  dff<-df
  
  dffx<-subset(df, select=-c(y,z))
  
  dfw<-df[ , !(names(df) %in% drop.x)]
  m1<-multinom(z~.,data = dfw[,-1],trace=FALSE)
  w<-numeric(nrow(dfw))
  pr<-predict(m1,type = "probs")
  #(table(dff$z)/n)[dff$z[i]]
  for (i in 1:nrow(dfw)) {
    w[i]<-1/pr[i,dff$z[i]]
  }
  dff$w<-w
  
  
  #Scores
  Sa<-predict(las1,s=las1$lambda.min,newx=as.matrix(dffx))-apply(cbind(predict(las2,s=las2$lambda.min,newx=as.matrix(dffx)),predict(las3,s=las3$lambda.min,newx=as.matrix(dffx))),1,max)
  Sb<-predict(las2,s=las2$lambda.min,newx=as.matrix(dffx))-apply(cbind(predict(las1,s=las1$lambda.min,newx=as.matrix(dffx)),predict(las3,s=las3$lambda.min,newx=as.matrix(dffx))),1,max)
  Sc<-predict(las3,s=las3$lambda.min,newx=as.matrix(dffx))-apply(cbind(predict(las2,s=las2$lambda.min,newx=as.matrix(dffx)),predict(las1,s=las1$lambda.min,newx=as.matrix(dffx))),1,max)
  
  dff$S<-apply(cbind(Sa,Sb,Sc),1,max)
  dff$d<-factor(apply(cbind(Sa,Sb,Sc),1,which.max),levels=c(1,2,3))
  
  ind1<-which(table(subset(dff,z==1)$d)>4)
  ind2<-which(table(subset(dff,z==2)$d)>4)
  ind3<-which(table(subset(dff,z==3)$d)>4)
  
  
  dff1<-subset(subset(dff,z==1),d%in%c(ind1))
  dff2<-subset(subset(dff,z==2),d%in%c(ind2))
  dff3<-subset(subset(dff,z==3),d%in%c(ind3))
  
  h1<-h.amise(dff1$S,deriv.order = 1)$h
  h2<-h.amise(dff2$S,deriv.order = 1)$h
  h3<-h.amise(dff3$S,deriv.order = 1)$h
  
  #Score model
  #my<-lm(y~S*d,data = dff)
  my1<-lm(y~S*d,weights = w,data = dff1)
  my2<-lm(y~S*d,weights = w,data = dff2)
  my3<-lm(y~S*d,weights = w,data = dff3)
  
  S0a<-predict(las1,s=las1$lambda.min,newx=as.matrix(xn))-apply(cbind(predict(las2,s=las2$lambda.min,newx=as.matrix(xn)),predict(las3,s=las3$lambda.min,newx=as.matrix(xn))),1,max)
  S0b<-predict(las2,s=las2$lambda.min,newx=as.matrix(xn))-apply(cbind(predict(las1,s=las1$lambda.min,newx=as.matrix(xn)),predict(las3,s=las3$lambda.min,newx=as.matrix(xn))),1,max)
  S0c<-predict(las3,s=las3$lambda.min,newx=as.matrix(xn))-apply(cbind(predict(las2,s=las2$lambda.min,newx=as.matrix(xn)),predict(las1,s=las1$lambda.min,newx=as.matrix(xn))),1,max)
  S0<-apply(cbind(S0a,S0b,S0c),1,max)
  d0<-factor(apply(cbind(S0a,S0b,S0c),1,which.max),levels = c(1,2,3))
  
  check<-data.frame(matrix(nrow=m,ncol = 3))
  tt<-data.frame(matrix(nrow=m,ncol = 2))
  nw.lam<-data.frame(matrix(nrow=m,ncol = 2))
  ll.lam<-data.frame(matrix(nrow=m,ncol = 2))
  for(k in 1:m){
    
    nonpar1<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=1,h1)
    nonpar2<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=2,h2)
    nonpar3<-alpha.hat.w(dff,si=S0[k],di=d0[k],zi=3,h3)
    
    ynw1<-nonpar1$nwe
    ynw2<-nonpar2$nwe
    ynw3<-nonpar3$nwe
    
    yll1<-nonpar1$lle
    yll2<-nonpar2$lle
    yll3<-nonpar3$lle
    
    ynw<-sort(c(ynw1,ynw2,ynw3))
    nw.lam[k,]<-c(ynw[3]-ynw[2],ynw[3]-ynw[1])
    
    yll<-sort(c(yll1,yll2,yll3))
    ll.lam[k,]<-c(yll[3]-yll[2],yll[3]-yll[1])
    
    check[k,]<-c(yn[k],which.max(c(ynw1,ynw2,ynw3)),which.max(c(yll1,yll2,yll3)))
    tt[k,1]<-ifelse(check[k,1]==check[k,2],1,0)
    tt[k,2]<-ifelse(check[k,1]==check[k,3],1,0)
  }
  C.nw<-max(nw.lam)
  C.ll<-max(ll.lam)
  
  nw.Lam12<-apply(nw.lam,2,mean)/C.nw
  ll.Lam12<-apply(ll.lam,2,mean)/C.ll
  nw.assign<-table(yn,check[,2])
  ll.assign<-table(yn,check[,3])
  accprob<-apply(tt,2,mean)
  return(list(nw.assign=nw.assign,ll.assign=ll.assign,accprob=accprob,nw.Lam12=nw.Lam12,ll.Lam12=ll.Lam12))
}


statin<-read.csv("D:\\UofL\\Project\\Data\\Chathura\\statin.csv") #The data set is available online in the Dryad Digital Repositary via:http://doi.org/10.5061/dryad.6qs53

statin$Ethnicity<-as.factor(statin$Ethnicity)
statin$Education<-as.factor(statin$Education)
statin$Albuminuria<-as.factor(statin$Albuminuria)
statin$z<-as.factor(statin$z)

statin<-na.omit(statin)


var_dum<-data.frame(Ethnicity=statin$Ethnicity,Education=statin$Education,Albuminuria=statin$Albuminuria)
results <- fastDummies::dummy_cols(var_dum, remove_first_dummy = TRUE)[,-c(1,2,3)]
dat<-data.frame(statin[,-c(4,5,15,18)],results,z=statin$z)
dat<-statin

#Preparing data for 5-fold cross validation
#Subset the data into 3 treatment groups
dat1<-subset(dat,z==1)
dat2<-subset(dat,z==2)
dat3<-subset(dat,z==3)


#FUntion to split a data set into 5 random sub groups
spec = c(rep(0.2,5))
datsp<-function(df){
  g = sample(cut(
    seq(nrow(df)), 
    nrow(df)*cumsum(c(0,spec)),
    labels = names(spec)
  ))
  res = split(df, g)
  names(res)<-c(1:5)
  return(res)
}

datsp1<-datsp(dat1)
datsp2<-datsp(dat2)
datsp3<-datsp(dat3)

train1<-rbind(bind_rows(datsp1[1:4]),bind_rows(datsp2[1:4]),bind_rows(datsp3[1:4]))
test1<-rbind(bind_rows(datsp1[5]),bind_rows(datsp2[5]),bind_rows(datsp3[5]))

train2<-rbind(bind_rows(datsp1[c(1,2,3,5)]),bind_rows(datsp2[c(1,2,3,5)]),bind_rows(datsp3[c(1,2,3,5)]))
test2<-rbind(bind_rows(datsp1[4]),bind_rows(datsp2[4]),bind_rows(datsp3[4]))

train3<-rbind(bind_rows(datsp1[c(1,2,4,5)]),bind_rows(datsp2[c(1,2,4,5)]),bind_rows(datsp3[c(1,2,4,5)]))
test3<-rbind(bind_rows(datsp1[3]),bind_rows(datsp2[3]),bind_rows(datsp3[3]))

train4<-rbind(bind_rows(datsp1[c(1,3,4,5)]),bind_rows(datsp2[c(1,3,4,5)]),bind_rows(datsp3[c(1,3,4,5)]))
test4<-rbind(bind_rows(datsp1[2]),bind_rows(datsp2[2]),bind_rows(datsp3[2]))

train5<-rbind(bind_rows(datsp1[c(2,3,4,5)]),bind_rows(datsp2[c(2,3,4,5)]),bind_rows(datsp3[c(2,3,4,5)]))
test5<-rbind(bind_rows(datsp1[1]),bind_rows(datsp2[1]),bind_rows(datsp3[1]))

#Function to obtain the optimal treatment assignment for a given testing data set using training data.
resultgbm<-function(test,train){
  m=nrow(test)
  
  df<-train
  xn<-test[,-c(1,ncol(test))]
  yn<-test$z
  return(match.gbm(df,xn))
}

resultgbm(test1,train1)
