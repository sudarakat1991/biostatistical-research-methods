library(robustHD)
library(survival)

xt.f<-function(t){
  a <- runif(1)
  b <- runif(1)
  w <- runif(1,0,2*pi)
  ret <- sin(w*t)+(a+2*pi)*t+b
  return(ret)
}

cgen <- function(n,c0){
  C <- runif(n,c0,14)
  return (C)
}

dataGen <- function(n, tj, c0){
  Xt <-array(NA,dim = c(n,length(tj)))
  yx<-numeric(n)
  for (i in 1:n) {
    Xt[i,]<-xt.f(tj)
    yx[i]<-2*mean(abs(Xt[i,])*(1-cos(pi*tj)))+rnorm(1,0,2)
  }
  C= cgen(n,c0)
  y<-apply(cbind(yx,C), 1,min)
  delt<-as.numeric(yx<C) 
  
  survdat.forsimulationonly = data.frame(yx=yx,C=C);
  survdat <-data.frame(y=y,delt=delt)
  survdat$Xt=Xt
  return(list(survdat.forsimulationonly, survdat))
}


surv.t<-function(data.train, xn, tvec, h){
  # data.train: observed data
  # xn:  the new X
  # tvec: the vector of time point of interest
  # h: the bandwidth
  
  #if (t>max(y)) { return(0) }
  y = data.train$y           # y: the observed time
  delt = data.train$delt     # delt: the censoring indicators
  Xtr = data.train$Xt        # Xtr: the observed X
  
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
        numer <-  Bx[j]/gx #1-((y-y[j])<=0)%*%(Bx)/gx
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

orig.surv<-function(tvec,x){
  # tvec: vector of time points
  # x: covariate
  surv.prob <-1-pnorm((tvec-2*mean(abs(x)*(1-cos(pi*tj)))),0,2)
  return(surv.prob)
}


### prediction error with true survival prob
prederr<-function(fitted, Xt, tvec){
  # fitted: fitted survival probability
  # Xt: new covariates where the true survival prob, must be a matrix
  # tvec: vector of time points
  # xn: the new covariate where the survival prob will be calculated
  
  # obtain true probabilities at tvec given a new covariate 
  suvt = 0*fitted
  for (j in 1:nrow(Xt)) {
    xn<-Xt[j,]
    suvt[j, ]<-  orig.surv(tvec,xn)
  }
  reserr = sqrt(apply((suvt - fitted)^2, 1, mean))
  return(reserr)
}


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


CVh <- function(df, K, tvec, h.grid){
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
      
      hmat.err.trueprob[,k] <- prederr(surv.prob.test, data.test$Xt, tvec)
      hmat.err.brier[,k] <- Brier.s(surv.prob.test, data.train, data.test, tvec, h)
    }
    hvec.err.trueprob[j] <- mean(hmat.err.trueprob,na.rm=T)
    hvec.err.brier[j] <- mean(hmat.err.brier,na.rm=T) #,na.rm =TRUE
  }
  #return(cbind(h.grid,hvec.err.trueprob,hvec.err.brier))
  return(list(hvec.err.trueprob,hvec.err.brier,c(h.grid[which.min(hvec.err.trueprob )], h.grid[which.min(hvec.err.brier)])))
}

##################################


nt = 100
tj<-seq(-1,1,by=0.01)
tvec <- seq(3,15,0.1)
h.grid<-seq(0.5,2,0.05)

#Testing optimal bandwidths: Censoring around 15%
result.n.c<-function(xx,n,c0){
  
  train.data1 <- dataGen(n, tj,c0)
  survdat1 <- train.data1[[2]]
  censor1<-mean(survdat1$delt)
  
  
  test.data1 <- dataGen(nt, tj,c0)
  survdat.test1 <- test.data1[[2]]
  
  tmpCV1=CVh(survdat1, K=5, tvec, h.grid)
  
  #Kaplan meier estimates for comparison
  fit<-survfit(Surv(survdat1$y, survdat1$delt) ~ 1)
  survest <- stepfun(fit$time, c(1, fit$surv))
  kmsurv<-survest(tvec)
  
  hk.err11<-array(NA,dim = c(nt,2))
  
  for (i in 1:nt) {
    xnt<-survdat.test1$Xt[i,]
    #h1<-8.45
    #h2<-0.48
    h<-tmpCV1[[3]][1]
    #h<-5.054688
    suvt<-surv.t(survdat1,xnt,tvec,h)
    samp<-orig.surv(tvec,xnt)
    hk.err11[i,]<-c(sqrt(mean((suvt - samp)^2)),sqrt(mean((kmsurv - samp)^2)))
  }
  
  hk.err12<-array(NA,dim = c(nt,2))
  
  for (i in 1:nt) {
    xnt<-survdat.test1$Xt[i,]
    #h1<-8.45
    #h2<-0.48
    h<-tmpCV1[[3]][2]
    #h<-5.054688
    suvt<-surv.t(survdat1,xnt,tvec,h)
    samp<-orig.surv(tvec,xnt)
    hk.err12[i,]<-c(sqrt(mean((suvt - samp)^2)),sqrt(mean((kmsurv - samp)^2)))
  }
  
  return(list(censor1,tmpCV1,hk.err11,hk.err12))
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

pars <- data.frame(xx=seq(1:8),n=rep(c(100,100,400,400),2),c0=rep(c(12,10,12,10),2))

start_time <- proc.time()
study_results <- clusterMap(my_cluster,result.n.c,xx=pars$xx,n=pars$n,c0=pars$c0)
print(proc.time() - start_time)


store.n1.c15<-study_results[[1]]
store.n1.c25<-study_results[[2]]
store.n2.c15<-study_results[[3]]
store.n2.c25<-study_results[[4]]


true10015<-c(store.n1.c15[[1]],(store.n1.c15[[2]][[3]])[1],apply(store.n1.c15[[3]],2,mean),apply(store.n1.c15[[3]],2,sd),store.n1.c15[[2]][[1]])
brier10015<-c(store.n1.c15[[1]],(store.n1.c15[[2]][[3]])[2],apply(store.n1.c15[[4]],2,mean),apply(store.n1.c15[[4]],2,sd),store.n1.c15[[2]][[2]])


true40015<-c(store.n2.c15[[1]],(store.n2.c15[[2]][[3]])[1],apply(store.n2.c15[[3]],2,mean),apply(store.n2.c15[[3]],2,sd),store.n2.c15[[2]][[1]])
brier40015<-c(store.n2.c15[[1]],(store.n2.c15[[2]][[3]])[2],apply(store.n2.c15[[4]],2,mean),apply(store.n2.c15[[4]],2,sd),store.n2.c15[[2]][[2]])

true10025<-c(store.n1.c25[[1]],(store.n1.c25[[2]][[3]])[1],apply(store.n1.c25[[3]],2,mean),apply(store.n1.c25[[3]],2,sd),store.n1.c25[[2]][[1]])
brier10025<-c(store.n1.c25[[1]],(store.n1.c25[[2]][[3]])[2],apply(store.n1.c25[[4]],2,mean),apply(store.n1.c25[[4]],2,sd),store.n1.c25[[2]][[2]])


true40025<-c(store.n2.c25[[1]],(store.n2.c25[[2]][[3]])[1],apply(store.n2.c25[[3]],2,mean),apply(store.n2.c25[[3]],2,sd),store.n2.c25[[2]][[1]])
brier40025<-c(store.n2.c25[[1]],(store.n2.c25[[2]][[3]])[2],apply(store.n2.c25[[4]],2,mean),apply(store.n2.c25[[4]],2,sd),store.n2.c25[[2]][[2]])


true1<-rbind(true10015,true10025,true40015,true40025)
brier1<-rbind(brier10015,brier10025,brier40015,brier40025)


#######################################################

store.n1.c15<-study_results[[5]]
store.n1.c25<-study_results[[6]]
store.n2.c15<-study_results[[7]]
store.n2.c25<-study_results[[8]]


true10015<-c(store.n1.c15[[1]],(store.n1.c15[[2]][[3]])[1],apply(store.n1.c15[[3]],2,mean),apply(store.n1.c15[[3]],2,sd),store.n1.c15[[2]][[1]])
brier10015<-c(store.n1.c15[[1]],(store.n1.c15[[2]][[3]])[2],apply(store.n1.c15[[4]],2,mean),apply(store.n1.c15[[4]],2,sd),store.n1.c15[[2]][[2]])


true40015<-c(store.n2.c15[[1]],(store.n2.c15[[2]][[3]])[1],apply(store.n2.c15[[3]],2,mean),apply(store.n2.c15[[3]],2,sd),store.n2.c15[[2]][[1]])
brier40015<-c(store.n2.c15[[1]],(store.n2.c15[[2]][[3]])[2],apply(store.n2.c15[[4]],2,mean),apply(store.n2.c15[[4]],2,sd),store.n2.c15[[2]][[2]])

true10025<-c(store.n1.c25[[1]],(store.n1.c25[[2]][[3]])[1],apply(store.n1.c25[[3]],2,mean),apply(store.n1.c25[[3]],2,sd),store.n1.c25[[2]][[1]])
brier10025<-c(store.n1.c25[[1]],(store.n1.c25[[2]][[3]])[2],apply(store.n1.c25[[4]],2,mean),apply(store.n1.c25[[4]],2,sd),store.n1.c25[[2]][[2]])


true40025<-c(store.n2.c25[[1]],(store.n2.c25[[2]][[3]])[1],apply(store.n2.c25[[3]],2,mean),apply(store.n2.c25[[3]],2,sd),store.n2.c25[[2]][[1]])
brier40025<-c(store.n2.c25[[1]],(store.n2.c25[[2]][[3]])[2],apply(store.n2.c25[[4]],2,mean),apply(store.n2.c25[[4]],2,sd),store.n2.c25[[2]][[2]])


true2<-rbind(true10015,true10025,true40015,true40025)
brier2<-rbind(brier10015,brier10025,brier40015,brier40025)

truef<-rbind(true1,true2)
brierf<-rbind(brier1,brier2)

write(t(truef),"true1.txt",ncolumns = ncol(truef))
write(t(brierf),"brier1.txt",ncolumns = ncol(brierf))