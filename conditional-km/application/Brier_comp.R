library(locpol)
library(MASS)
library(survival)
library(survAUC)
library(survcomp)

setwd("D:\\UofL\\Project\\Dr.Qi\\Read\\Cox Regression\\FLCRM\\biom12748-sup-0002-suppinfo-s2")

load("hippocampussurface.dat")
load("clinical.dat")

#block<-as.factor(read.table("../block.txt")[,1])
#h.grid1<-read.table("sud.txt")[,]

fulldf<-data.frame(y=scoredata[,1],delt=scoredata[,2])
#fulldf$Xt<-t(cleandata)

xmean<-mean(scoredata[,15])
xsd<-sd(scoredata[,15])

fulldf$X<-scale(scoredata[,15],center = xmean, scale=xsd)

xmeanf<-apply(t(cleandata),2,mean)
xsdf<-rep(sd(apply(t(cleandata),1,function (x) sqrt(sum(x^2)))),length(xmeanf))
#xsd<-apply(Xtr,2,sd)

fulldf$Xt<-scale(t(cleandata),center =xmeanf,scale = xsdf)

#set.seed(1)
ind<-seq(1:273)
traindf<-fulldf[ind,]
testdf<-fulldf[-ind,]

tvec <- seq(100,1600,10)
n<-nrow(traindf)
nt<-nrow(testdf)
hf<-0.01


scaleX<-scale(traindf$Xt,center=TRUE,scale=FALSE)
Eigvec=svd(scaleX)$v
Singval=svd(scaleX)$d
est_Eigval=Singval^2*hf/(n-1)

totaleigenscore<-scaleX%*%(Eigvec/sqrt(hf))*hf

ordertime<-order(traindf$y)
sorttime<-traindf$y[ordertime]
sortevent<-traindf$delt[ordertime]

sortfailureindex<-which(sortevent==1)

hh<-40

#### AIC selection result
AIC<-loglikelihood<-rep(0,hh)
for (ss in 1:hh) {
  est_Eigfun=Eigvec[,1:ss]/sqrt(hf)
  eigenscore<-scaleX%*%est_Eigfun*hf
  designmatrix<-cbind(eigenscore,traindf$X)
  fullcoxresult<-coxph(formula=Surv(traindf$y, traindf$delt)~designmatrix)
  orderdesignmatrix<-designmatrix[ordertime,]
  estimatecoefficient<-fullcoxresult$coefficients
  fullcoxcurve<-survfit(fullcoxresult)
  finalresult<-summary(fullcoxresult)
  loglikelihood[ss]<-finalresult$loglik[2]
  AIC[ss]<-2*ss-2*loglikelihood[ss]
}

plot(AIC)

AICselectindex<-which.min(AIC)

AICselectindex<-18

est_Eigfun=Eigvec[,1:AICselectindex]/sqrt(hf)
eigenscore<-scaleX%*%est_Eigfun*hf
designmatrix<-cbind(eigenscore,traindf$X)

newdataframe<-data.frame(time=traindf$y, event=traindf$delt, designmatrix)
fullcoxresult<-coxph(formula=Surv(time,event)~., data=newdataframe)

finalresult<-summary(fullcoxresult)
summarycoefficient<-finalresult$coefficients
coefficient<-summarycoefficient[,1]
secoefficient<-summarycoefficient[,3]

#### Estimated Beta
AICestimatebeta<-est_Eigfun%*%as.matrix(coefficient[1:AICselectindex])

#### Estimated Gamma
AICestimategamma<-coefficient[-(1:AICselectindex)]

#Estimate
scaleX.t<-scale(testdf$Xt,center=TRUE,scale=FALSE)
Eigvec.t=svd(scaleX.t)$v
Singval.t=svd(scaleX.t)$d
est_Eigval.t=Singval.t^2*hf/(nt-1)

est_Eigfun.t=Eigvec.t[,1:AICselectindex]/sqrt(hf)
eigenscore.t<-scaleX.t%*%est_Eigfun.t*hf
designmatrix.t<-cbind(eigenscore.t,testdf$X)
colnames(designmatrix.t)<-colnames(newdataframe)[-c(1,2)]


fullcoxcurve<-survfit(fullcoxresult,newdata = designmatrix.t[1,])
survest.cox <- stepfun(fullcoxcurve$time, c(1, fullcoxcurve$surv))
coxest<-survest.cox(tvec)
#################################################

#####Censored
data.C <- traindf
data.C$delt <- 1-data.C$delt

sortevent.c<-data.C$delt[ordertime]
sortfailureindex.c<-which(sortevent.c==1)

newdataframe.c<-data.frame(time=traindf$y, event=data.C$delt, designmatrix)
fullcoxresult.c<-coxph(formula=Surv(time,event)~., data=newdataframe.c)

surv.t<-function(data.train, xnt, xn, tvec, h1,h2){
  # data.train: observed data
  # xn:  the new X
  # tvec: the vector of time point of interest
  # h: the bandwidth
  
  #if (t>max(y)) { return(0) }
  y = data.train$y           # y: the observed time
  delt = data.train$delt     # delt: the censoring indicators
  Xtr = data.train$Xt        # Xtr: the observed X
  X = data.train$X
  n <- length(y)
  
  
  difxnorm1 <- apply(scale(Xtr,center =xnt,scale = FALSE ), 
                     1, function(x) sqrt(sum(x^2)))
  Bx1 <- dnorm(difxnorm1/h1)/(n*h1)
  Bx2 <- dnorm((X-xn)/h2)/(h2)
  Bx <- Bx1*Bx2
  
  # gz in Dabrowska 1989, the density of X at xn
  #gx1<-sum(Bx1)              
  #gx2<-sum(Bx2)
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
        numer <- (Bx[j]/gx)
        denom <- 1-((y-y[j])<0)%*%((Bx)/gx)
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

h1=1.1;h2=0.6

pit.FLCRM <- matrix(0, nrow(testdf), length(tvec))
pit.FKM <- matrix(0, nrow(testdf), length(tvec))
pit.KM <- matrix(0, nrow(testdf), length(tvec))

for (i in 1:nt) {
  fullcoxcurve<-survfit(fullcoxresult,newdata = designmatrix.t[i,])
  survest.cox <- stepfun(fullcoxcurve$time, c(1, fullcoxcurve$surv))
  
  fit<-survfit(Surv(traindf$y, traindf$delt) ~ 1)
  survest <- stepfun(fit$time, c(1, fit$surv))
  
  pit.KM[i, ]<-survest(tvec)
  pit.FLCRM[i, ]<-survest.cox(tvec)
  pit.FKM[i, ] <- surv.t(traindf, testdf$Xt[i,], testdf$X[i], tvec, h1,h2)
}


fit.c<-survfit(Surv(traindf$y, data.C$delt) ~ 1)
survest.c <- stepfun(fit.c$time, c(1, fit.c$surv))


data.C <- traindf
data.C$delt <- 1-data.C$delt
  
nnew <- nrow(testdf)
wit.FKM <- matrix(0,nnew, length(tvec))
n.tilde.FKM<-0*wit.FKM

wit.KM <- matrix(0,nnew, length(tvec))
n.tilde.KM<-0*wit.KM

wit.FLCRM <- matrix(0,nnew, length(tvec))
n.tilde.FLCRM<-0*wit.FLCRM


for(i in 1:nnew){
    # for each row in data.new
    xnt <- testdf$Xt[i,]
    xn <- testdf$X[i]
    yn <- testdf$y[i]
    delt <- testdf$delt[i]
    # calculate the probability of P(C>yn|xn) and P(C>tvec|xn)
    
    fullcoxcurve.c<-survfit(fullcoxresult.c,newdata = designmatrix.t[i,])
    survest.cox.c <- stepfun(fullcoxcurve.c$time, c(1, fullcoxcurve.c$surv))
    
    probCy.FKM <-max(surv.t(data.C, xnt, xn, yn, h1, h2),1E-6)
    probCt.FKM <-apply(cbind(surv.t(data.C, xnt, xn, tvec, h1, h2),1E-6), 1,max)
    
    probCy.FLCRM <-max(survest.cox(tvec),1E-6)
    probCt.FLCRM <-apply(cbind(survest.cox(tvec),1E-6), 1,max)
    
    probCy.KM <-max(survest.c(tvec),1E-6)
    probCt.KM <-apply(cbind(survest.c(tvec),1E-6), 1,max)
    
    
    Itd<-1*(yn<=tvec)*(delt==1)
    yit<- 1*(yn>tvec)
    
    n.tilde.FKM[i,]<-Itd/probCy.FKM+yit/probCt.FKM
    wit.FKM[i,]<-((pit.FKM[i,]^2)*Itd)/probCy.FKM+(((1-pit.FKM[i,])^2)*yit)/probCt.FKM
    
    n.tilde.KM[i,]<-Itd/probCy.KM+yit/probCt.KM
    wit.KM[i,]<-((pit.KM[i,]^2)*Itd)/probCy.KM+(((1-pit.KM[i,])^2)*yit)/probCt.KM
    
    n.tilde.FLCRM[i,]<-Itd/probCy.FLCRM+yit/probCt.FLCRM
    wit.FLCRM[i,]<-((pit.FLCRM[i,]^2)*Itd)/probCy.FLCRM+(((1-pit.FLCRM[i,])^2)*yit)/probCt.FLCRM
    
}
Br.FKM<-apply(wit.FKM,2,sum,na.rm=T)/apply(n.tilde.FKM,2,sum,na.rm=T)
Br.KM<-apply(wit.KM,2,sum,na.rm=T)/apply(n.tilde.KM,2,sum,na.rm=T)
Br.FLCRM<-apply(wit.FLCRM,2,sum,na.rm=T)/apply(n.tilde.FLCRM,2,sum,na.rm=T)

Br.dat<-data.frame(Time=tvec,KM=Br.KM,Cond.KM=Br.FKM,FLCRM=Br.FLCRM)

Br.dat1<-Br.dat%>% gather(Method, Brier.Score, 2:4)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 3
mycols = gg_color_hue(n)

Br.dat1$Method<-as.factor(Br.dat1$Method)
levels(Br.dat1$Method)<-c("Cond.KM.Brier (AUC=219.5)", "FLCRM (AUC=326.7)", "KM (AUC=418.7)" )

figure0 <- ggplot(data = Br.dat1) + 
  geom_line(aes(x = Time , y = Brier.Score,col=Method),size=0.8)+
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12 ),
         axis.title = element_text( size = 14),
         legend.title = element_text("Method", size = 16 ),
         legend.text=element_text( size = 14 ),
         legend.position="bottom");figure0


library(DescTools)

AUC(x=Br.dat$Cond.KM, y=Br.dat$Time, method="trapezoid")
AUC(x=Br.dat$KM, y=Br.dat$Time, method="trapezoid")
AUC(x=Br.dat$FLCRM, y=Br.dat$Time, method="trapezoid")

sum(is.na(Br.dat))
