library(survival)
library(robustHD)
library(MASS)

### Generate data
xt.f<-function(t){
  a <- rnorm(1)
  b <- rnorm(1)
  v1 = rep(0, 10)
  v2 = rep(0, 10)
  tmp = rep(0,length(t))
  for(k in 1:10){
    v1[k] = rnorm(1,0,1/k)
    v2[k] = rnorm(1,0,1/k)   
    tmp = tmp + v1[k]*sin(2*(2*k-1)*pi*t) + v2[k]*cos(2*(2*k-1)*pi*t)
  }
  ret <- a + b*t + tmp
  return(ret)
}

cgen <- function(n,c0){
  # c0 is chosen to achieve a desired censoring rate of 10%, 30% or 50%
  C <- runif(n,0,c0)
  return (C)
}

dataGen <- function(n, tj, lam, c0){
  Xt <-array(NA,dim = c(n,length(tj)))
  yx<-numeric(n)
  for (i in 1:n) {
    Xt[i,]<-xt.f(tj)
    beta <-  0.6*(sin(pi*tj)-cos(pi*tj) +
                    sin(3*pi*tj/10) - cos(3*pi*tj) +
                    sin(5*pi*tj)/9 - cos(5*pi*tj)/9+
                    sin(7*pi*tj)/16 - cos(7*pi*tj)/16+
                    sin(9*pi*tj)/25 - cos(9*pi*tj)/25 + 
                    (2*pi)^(-1/2)* exp(-(tj-0.5)^2/2)) 
    U<-runif(1)
    yx[i]<--(log(U))/(lam*exp(mean(Xt[i,]*beta)))
  }
  C= cgen(n,c0)
  y<-apply(cbind(yx,C), 1,min)
  delt<-as.numeric(yx<C) 
  
  survdat.forsimulationonly = data.frame(yx=yx,C=C);
  survdat <-data.frame(y=y,delt=delt)
  survdat$Xt=Xt
  return(list(survdat.forsimulationonly, survdat))
}


### true survival prob
orig.surv<-function(tvec,x,beta, lam){
  # tvec: vector of time points
  # x: covariate
  surv.prob <-exp(-lam*tvec*exp(mean(beta*x)))
  return(surv.prob)
}



#Estimation: Survival func - KM for functional data
surv.t<-function(data.train, xn, tvec, h){
  # data.train: observed data
  # xn:  the new X
  # tvec: the vector of time point of interest
  # h: the bandwidth
  
  #if (t>max(y)) { return(0) }
  y = data.train$y           # y: the observed time
  delt = data.train$delt     # delt: the censoring indicators
  Xtr = data.train$Xt        # Xtr: the observed X
  
  #Xtr<-t(standardize(t(Xtr)))
  #xn<-standardize(xn)
  
  n <- length(y)
  
  xmean<-apply(Xtr,2,mean)
  xsd<-rep(sd(apply(Xtr,1,function (x) sqrt(sum(x^2)))),length(xmean))
  #xsd<-apply(Xtr,2,sd)
  
  Xtr<-scale(Xtr,center =xmean,scale = xsd)
  xn<-(xn-xmean)/xsd
  
  difxnorm <- apply(scale(Xtr,center =xn,scale = FALSE ), 
                    1, function(x) sqrt(sum(x^2)))
  
  Bx <- dnorm(difxnorm/h)/(n*h)
  
  # gz in Dabrowska 1989, the density of X at xn
  gx<-max(sum(Bx),1E-6)             
  
  # vector to store the estimated survival prob at tvec given xn
  surv.prob.est <- rep(0, length(tvec))
  
  for(i in 1:length(tvec)){
    # for all time points in the vector tvec
    t = tvec[i]
    
    fac <- as.numeric(length(y))
    # betat in W. Gonzalez-Manteiga & C. Cadarso-Suarez (1994)
    betat <- 1*(y<=t)*(delt==1)
    
    for (j in 1:length(y)) {
      I <- betat[j]
      if (I==0) {
        fac[j]=1
      } else {
        numer <- Bx[j]/gx
        denom <- 1-((y-y[j])<0)%*%(Bx)/gx
        if (denom==0) {
          fac[j]=0
        } else {
          fac[j]=exp(-numer/denom)
        }
      }
    }
    surv.prob.est[i]<- prod(fac)
  }
  
  return(surv.prob.est)
}


### prediction error with true survival prob
prederr<-function(fitted, Xt, tvec,lam){
  # fitted: fitted survival probability
  # Xt: new covariates where the true survival prob, must be a matrix
  # tvec: vector of time points
  # xn: the new covariate where the survival prob will be calculated
  
  # obtain true probabilities at tvec given a new covariate 
  suvt = 0*fitted
  for (j in 1:nrow(Xt)) {
    suvt[j, ]<-  orig.surv(tvec,Xt[j,],beta,lam)
  }
  reserr = sqrt(apply((suvt - fitted)^2, 1, mean))
  return(reserr)
}


#Brier Score
Brier.s<-function(pit,data.train,data.new,tvec,h){
  data.C <- data.train
  data.C$delt <- 1-data.C$delt
  
  nnew <- nrow(data.new)
  wit <- matrix(0,nnew, length(tvec))
  n.tilde<-0*wit
  for(i in 1:nnew){
    # for each row in data.new
    xn = data.new$Xt[i,]
    yn = data.new$y[i]
    delt = data.new$delt[i]
    # calculate the probability of P(C>yn|xn) and P(C>tvec|xn)
    
    probCy <-max(surv.t(data.C, xn, yn, h),1E-6)
    probCt <-apply(cbind(surv.t(data.C, xn, tvec, h),1E-6), 1,max)
    
    Itd<-1*(yn<=tvec)*(delt==1)
    yit<- 1*(yn>tvec)
    n.tilde[i,]<-Itd/probCy+yit/probCt
    wit[i,]<-((pit[i,]^2)*Itd)/probCy+(((1-pit[i,])^2)*yit)/probCt
    
  }
  Br<-apply(wit,2,sum,na.rm=T)/apply(n.tilde,2,sum,na.rm=T)
  return(Br)
}



### Cross-Validation
# K fold cross-validation
CVh <- function(df, K, tvec, h.grid,lam){
  # df: the data frame
  # K: K fold
  # tvec: time points of interest
  # h.grid: the grid for the bandwidth
  n = nrow(df)
  
  # get K folds
  # size: the size of testing set
  size <- n %/% K; 
  rdm <- runif(n);
  ranked <- rank(rdm);
  block <- (ranked-1) %/% size+1;
  block <- as.factor(block)
  
  # vectors to store cv errors using the true probabilities and brier scores
  hvec.err.trueprob <- rep(0,length(h.grid))
  hvec.err.brier <- rep(0,length(h.grid))
  
  for (j in 1:length(h.grid)) {
    h<-h.grid[j]
    ### matrix to store the cv errors using the true probabilities and brier scores for each fold
    hmat.err.trueprob <- array(NA,dim = c(size,K))
    hmat.err.brier <-  array(NA,dim = c(length(tvec),K))
    
    for(k in 1:K){
      ### training set and testing set
      data.train = df[which(block!=k),]
      data.test = df[which(block==k),]
      
      ### estimated survival probability for testing set
      surv.prob.test <- matrix(0, nrow(data.test), length(tvec))
      for (i in 1:size) {
        surv.prob.test[i, ] <- surv.t(data.train, data.test$Xt[i,], tvec, h)
      }
      
      hmat.err.trueprob[,k] <- prederr(surv.prob.test, data.test$Xt, tvec,lam)
      hmat.err.brier[,k] <- Brier.s(surv.prob.test, data.train, data.test, tvec, h)
    }
    
    hvec.err.trueprob[j] <- mean(hmat.err.trueprob,na.rm=T)
    hvec.err.brier[j] <- mean(hmat.err.brier,na.rm=T)
  }
  #return(cbind(h.grid,hvec.err.trueprob,hvec.err.brier))
  return(list(hvec.err.trueprob,hvec.err.brier,c(h.grid[which.min(hvec.err.trueprob )], h.grid[which.min(hvec.err.brier)])))
}


nt = 100
tj<-seq(0,1,by=0.01)
tvec <- seq(0.05,8,0.05)
h.grid<-seq(0.5,2,by=0.05)

beta <-  0.6*(sin(pi*tj)-cos(pi*tj) +
                sin(3*pi*tj/10) - cos(3*pi*tj) +
                sin(5*pi*tj)/9 - cos(5*pi*tj)/9+
                sin(7*pi*tj)/16 - cos(7*pi*tj)/16+
                sin(9*pi*tj)/25 - cos(9*pi*tj)/25 + 
                (2*pi)^(-1/2)* exp(-(tj-0.5)^2/2)) 

result.n.c<-function(xx,n,c0){
  train.data5 <- dataGen(n, tj,lam=1,c0)
  survdat5 <- train.data5[[2]]
  censor5<-mean(survdat5$delt)
  
  test.data5 <- dataGen(nt, tj,lam = 1,c0)
  survdat.test5 <- test.data5[[2]]
  
  
  tmpCV5 = CVh(survdat5, K=5, tvec, h.grid,lam = 1)
  
  ####################FLCRM
  hf<-tj[2]
  
  
  scaleX<-scale(survdat5$Xt,center=TRUE,scale=FALSE)
  Eigvec=svd(scaleX)$v
  Singval=svd(scaleX)$d
  est_Eigval=Singval^2*hf/(n-1)
  
  totaleigenscore<-scaleX%*%(Eigvec/sqrt(hf))*hf
  
  ordertime<-order(survdat5$y)
  sorttime<-survdat5$y[ordertime]
  sortevent<-survdat5$delt[ordertime]
  
  sortfailureindex<-which(sortevent==1)
  
  hh<-8
  
  #### AIC selection result
  # AIC<-loglikelihood<-rep(0,hh)
  # for (ss in 1:hh) {
  #   est_Eigfun=Eigvec[,1:ss]/sqrt(hf)
  #   eigenscore<-scaleX%*%est_Eigfun*hf
  #   designmatrix<-eigenscore
  #   fullcoxresult<-coxph(formula=Surv(survdat5$y, survdat5$delt)~designmatrix)
  #   orderdesignmatrix<-designmatrix[ordertime,]
  #   estimatecoefficient<-fullcoxresult$coefficients
  #   fullcoxcurve<-survfit(fullcoxresult)
  #   finalresult<-summary(fullcoxresult)
  #   loglikelihood[ss]<-finalresult$loglik[2]
  #   AIC[ss]<-2*ss-2*loglikelihood[ss]
  # }
  
  
  #AICselectindex<-which.min(AIC)
  AICselectindex<-10
  
  est_Eigfun=Eigvec[,1:AICselectindex]/sqrt(hf)
  eigenscore<-scaleX%*%est_Eigfun*hf
  designmatrix<-eigenscore
  
  newdataframe<-data.frame(time=survdat5$y, event=survdat5$delt, designmatrix)
  fullcoxresult<-coxph(formula=Surv(time,event)~., data=newdataframe)
  
  finalresult<-summary(fullcoxresult)
  summarycoefficient<-finalresult$coefficients
  coefficient<-summarycoefficient[,1]
  secoefficient<-summarycoefficient[,3]
  
  #### Estimated Beta
  AICestimatebeta<-est_Eigfun%*%as.matrix(coefficient[1:AICselectindex])
  
  scaleX.t<-scale(survdat.test5$Xt,center=TRUE,scale=FALSE)
  Eigvec.t=svd(scaleX.t)$v
  Singval.t=svd(scaleX.t)$d
  est_Eigval.t=Singval.t^2*hf/(nt-1)
  
  est_Eigfun.t=Eigvec.t[,1:AICselectindex]/sqrt(hf)
  eigenscore.t<-scaleX.t%*%est_Eigfun.t*hf
  designmatrix.t<-eigenscore.t
  colnames(designmatrix.t)<-colnames(newdataframe)[-c(1,2)]
  
  
  #Kaplan meier estimates for comparison
  fit<-survfit(Surv(survdat5$y, survdat5$delt) ~ 1)
  survest <- stepfun(fit$time, c(1, fit$surv))
  kmsurv<-survest(tvec)
  
  
  #Testing optimal bandwidths
  hk.err51<-array(NA,dim = c(nt,4))
  for (i in 1:nt) {
    xnt<-survdat.test5$Xt[i,]
    h1<-(tmpCV5[[3]])[1] 
    h2<-(tmpCV5[[3]])[2]
    #h<-5.054688
    fullcoxcurve<-survfit(fullcoxresult,newdata = designmatrix.t[i,])
    survest.cox <- stepfun(fullcoxcurve$time, c(1, fullcoxcurve$surv))
    coxest<-survest.cox(tvec)
    suvt.tr<-surv.t(survdat5,xnt,tvec,h1)
    suvt.br<-surv.t(survdat5,xnt,tvec,h2)
    #coxbet<-orig.surv(tvec,xnt,AICestimatebeta, lam=1)
    samp<-orig.surv(tvec,xnt,beta, lam=1)
    hk.err51[i,]<-c(sqrt(mean((suvt.tr - samp)^2)),sqrt(mean((suvt.br - samp)^2)),sqrt(mean((kmsurv - samp)^2)),sqrt(mean((coxest - samp)^2)))
  }
  mean.MSE<-apply(hk.err51,2,mean)
  sd.MSE<-apply(hk.err51,2,sd)
  return(list(censor5,tmpCV5,mean.MSE,sd.MSE))
}


library(parallel)
num_cores <- 10 #must be 20 or less
my_cluster <- makeCluster(num_cores)

clusterEvalQ(my_cluster, {
  library(survival)
  library(robustHD)
  library(MASS)
})

clusterExport(my_cluster,varlist = ls())

pars <- data.frame(xx=seq(1:8),n=rep(c(100,100,400,400),2),c0=rep(c(10,5,10,5),2))

start_time <- proc.time()
study_results <- clusterMap(my_cluster,result.n.c,xx=pars$xx,n=pars$n,c0=pars$c0)
print(proc.time() - start_time)

store.n1.c15<-study_results[[1]]
store.n1.c25<-study_results[[2]]
store.n2.c15<-study_results[[3]]
store.n2.c25<-study_results[[4]]


true10015<-c(store.n1.c15[[1]],store.n1.c15[[2]][[3]][1],store.n1.c15[[3]][c(1,3,4)],store.n1.c15[[4]][c(1,3,4)],store.n1.c15[[2]][[1]])
brier10015<-c(store.n1.c15[[1]],store.n1.c15[[2]][[3]][2],store.n1.c15[[3]][c(2,3,4)],store.n1.c15[[4]][c(2,3,4)],store.n1.c15[[2]][[2]])

true40015<-c(store.n2.c15[[1]],store.n2.c15[[2]][[3]][1],store.n2.c15[[3]][c(1,3,4)],store.n2.c15[[4]][c(1,3,4)],store.n2.c15[[2]][[1]])  
brier40015<-c(store.n2.c15[[1]],store.n2.c15[[2]][[3]][2],store.n2.c15[[3]][c(2,3,4)],store.n2.c15[[4]][c(2,3,4)],store.n2.c15[[2]][[2]])

true10025<-c(store.n1.c25[[1]],store.n1.c25[[2]][[3]][1],store.n1.c25[[3]][c(1,3,4)],store.n1.c25[[4]][c(1,3,4)],store.n1.c25[[2]][[1]])
brier10025<-c(store.n1.c25[[1]],store.n1.c25[[2]][[3]][2],store.n1.c25[[3]][c(2,3,4)],store.n1.c25[[4]][c(2,3,4)],store.n1.c25[[2]][[2]])

true40025<-c(store.n2.c25[[1]],store.n2.c25[[2]][[3]][1],store.n2.c25[[3]][c(1,3,4)],store.n2.c25[[4]][c(1,3,4)],store.n2.c25[[2]][[1]])
brier40025<-c(store.n2.c25[[1]],store.n2.c25[[2]][[3]][2],store.n2.c25[[3]][c(2,3,4)],store.n2.c25[[4]][c(2,3,4)],store.n2.c25[[2]][[2]])

true1<-rbind(true10015,true10025,true40015,true40025)
brier1<-rbind(brier10015,brier10025,brier40015,brier40025)


#######################
store.n1.c15<-study_results[[5]]
store.n1.c25<-study_results[[6]]
store.n2.c15<-study_results[[7]]
store.n2.c25<-study_results[[8]]

true10015<-c(store.n1.c15[[1]],store.n1.c15[[2]][[3]][1],store.n1.c15[[3]][c(1,3,4)],store.n1.c15[[4]][c(1,3,4)],store.n1.c15[[2]][[1]])
brier10015<-c(store.n1.c15[[1]],store.n1.c15[[2]][[3]][2],store.n1.c15[[3]][c(2,3,4)],store.n1.c15[[4]][c(2,3,4)],store.n1.c15[[2]][[2]])

true40015<-c(store.n2.c15[[1]],store.n2.c15[[2]][[3]][1],store.n2.c15[[3]][c(1,3,4)],store.n2.c15[[4]][c(1,3,4)],store.n2.c15[[2]][[1]])  
brier40015<-c(store.n2.c15[[1]],store.n2.c15[[2]][[3]][2],store.n2.c15[[3]][c(2,3,4)],store.n2.c15[[4]][c(2,3,4)],store.n2.c15[[2]][[2]])

true10025<-c(store.n1.c25[[1]],store.n1.c25[[2]][[3]][1],store.n1.c25[[3]][c(1,3,4)],store.n1.c25[[4]][c(1,3,4)],store.n1.c25[[2]][[1]])
brier10025<-c(store.n1.c25[[1]],store.n1.c25[[2]][[3]][2],store.n1.c25[[3]][c(2,3,4)],store.n1.c25[[4]][c(2,3,4)],store.n1.c25[[2]][[2]])

true40025<-c(store.n2.c25[[1]],store.n2.c25[[2]][[3]][1],store.n2.c25[[3]][c(1,3,4)],store.n2.c25[[4]][c(1,3,4)],store.n2.c25[[2]][[1]])
brier40025<-c(store.n2.c25[[1]],store.n2.c25[[2]][[3]][2],store.n2.c25[[3]][c(2,3,4)],store.n2.c25[[4]][c(2,3,4)],store.n2.c25[[2]][[2]])

true2<-rbind(true10015,true10025,true40015,true40025)
brier2<-rbind(brier10015,brier10025,brier40015,brier40025)

truef<-rbind(true1,true2)
brierf<-rbind(brier1,brier2)

write(t(truef),"true1.txt",ncolumns = ncol(truef))
write(t(brierf),"brier1.txt",ncolumns = ncol(brierf))