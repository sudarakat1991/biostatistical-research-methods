library(ggplot2)
library(tidyr)
library(cowplot)
library(plyr)
library(ggpubr)
library(gridExtra)
###########Scenario1

setwd("D:\\UofL\\Project\\Dr.Qi\\Simulations\\Results\\Codes\\CV\\Simu\\standardized_newbr\\Scenario1\\")

h.grid0<-seq(0.5,2,0.05)
hn<-length(h.grid0)

true1<-read.table("true1.txt",na.strings = "NA")
brier1<-read.table("brier1.txt",na.strings = "NA")

colnames(true1)<-c("Censoring_rate","h.opt","cond.MSE.mean","KM.MSE.mean","Cond.MSE.sd","KM.MSE.sd",paste(rep("h",hn),1:hn,rep(".err",hn),sep = ""))
colnames(brier1)<-colnames(true1)

true1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
true1$cens<-rep(c("15% Censoring","25% Censoring"),200)

brier1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
brier1$cens<-rep(c("15% Censoring","25% Censoring"),200)

sumdat <- data.frame(samp=true1$samp,cens=true1$cens,Cond.KM.MSE=true1$cond.MSE.mean,Cond.KM.Brier=brier1$cond.MSE.mean,KM=brier1$KM.MSE.mean) %>% gather(Method, MSE, 3:5)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 4
mycols = gg_color_hue(n)

ph0<-ggplot(data = sumdat, aes(y=MSE, x=Method)) + 
  geom_boxplot(aes(fill=Method))+
  #scale_fill_manual(values=c("#AD0100", "#7cAE00", "#00BFC4")) +
  facet_wrap( ~ cens+samp,ncol = 2)+
  labs(y="MSE",x="Method",title = "Scenario 1")+
  theme(legend.position = "bottom",axis.title.x=element_blank())+
  scale_fill_manual(values=mycols[-4])+
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12,angle = 15, hjust=1 ),
         axis.title = element_text( size = 16),
         legend.title = element_text( size = 16 ),
         legend.text=element_text( size = 14 ),
         strip.text = element_text(size = 16),
         plot.title = element_text(size = 18,hjust = 0.5));ph0

###############Chosing the optimal bandwidth

true1<-subset(true1,samp=="n=100")
brier1<-subset(brier1,samp=="n=100")

h.true0<-data.frame(h.grid0,t(aggregate(true1[,(ncol(true1)-hn-2+1):(ncol(true1)-2)],by=list(factor(true1$samp),factor(true1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))
h.brier0<-data.frame(h.grid0,t(aggregate(brier1[,(ncol(brier1)-hn-2+1):(ncol(brier1)-2)],by=list(factor(brier1$samp),factor(brier1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))

colnames(h.true0)<-c("h.grid","n100_c15","n100_c25")
colnames(h.brier0)<-c("h.grid","n100_c15","n100_c25")

h.true0<-h.true0%>% gather(Method, Error, 2:3)
h.brier0<-h.brier0%>% gather(Method, Error, 2:3)

opt.htrue<-data.frame(samp=factor(true1$samp),cens=factor(true1$cens),h.opt=true1$h.opt)
vert.htrue<-data.frame(Method=levels(factor(h.brier0$Method)),h.opt=ddply(opt.htrue,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

opt.hbrier<-data.frame(samp=factor(brier1$samp),cens=factor(brier1$cens),h.opt=brier1$h.opt)
vert.hbrier<-data.frame(Method=levels(factor(h.brier0$Method)),h.opt=ddply(opt.hbrier,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

face_lab<-as_labeller(c("n100_c15"="n=100, CR=15%", "n100_c25"="n=100, CR=25%"))


figure0 <- ggplot(data = h.brier0) + 
  geom_line(aes(x = h.grid , y = Error))+
  geom_vline(data = vert.htrue, aes(xintercept=h.opt,linetype="Brier"),size=0.8)+
  geom_vline(data = vert.hbrier, aes(xintercept=h.opt,linetype="MSE"),size=0.8)+
  #ylim(0.092,0.098)+
  geom_text(x=3, y=30, label="Scatter plot")+
  facet_wrap( ~ Method,ncol=1,labeller = face_lab)+
  labs(y="Brier Score",x="h.grid",title = "Scenario 1")+
  #scale_y_continuous(breaks=seq(0.07,0.15,by=0.0025))+
  theme(legend.position = "bottom")+
  #geom_text(data = ann_text, mapping= aes(x = 1.5, y=0.0975,label = label),size=5)+
  scale_color_manual(name = "Bandwidth", values = c(MSE = "#F8766D", Brier = "#7CAE00"))+
  scale_linetype_manual(name = "Bandwidth", values = c(MSE = "dashed", Brier = "dotted"))+
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12 ),
         axis.title = element_text( size = 14),
         legend.title = element_text( size = 14 ),
         legend.text=element_text( size = 12 ),
         strip.text = element_text(size = 14),
         plot.title = element_text(size = 14,hjust = 0.5));figure0


###########Scenario 2

setwd("D:\\UofL\\Project\\Dr.Qi\\Simulations\\Results\\Codes\\CV\\Simu\\standardized_newbr\\Scenario2_re\\")

h.grid11<-seq(0.5,2,by=0.1)
h.grid12<-seq(0.2,1,by=0.05)

hn<-length(h.grid11)

true1<-read.table("true1.txt",na.strings = "NA")
brier1<-read.table("brier1.txt",na.strings = "NA")


colnames(true1)<-c("Censoring_rate","h.opt","h2.opt","cond.MSE.mean","KM.MSE.mean","Cond.MSE.sd","KM.MSE.sd",paste(rep("h",hn),1:hn,rep(".err",hn),sep = ""))
colnames(brier1)<-colnames(true1)

true1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
true1$cens<-rep(c("25% Censoring","15% Censoring"),200)

brier1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
brier1$cens<-rep(c("25% Censoring","15% Censoring"),200)

sumdat1 <- data.frame(samp=true1$samp,cens=true1$cens,Cond.KM.MSE=true1$cond.MSE.mean,Cond.KM.Brier=brier1$cond.MSE.mean,KM=brier1$KM.MSE.mean) %>% gather(Method, MSE, 3:5)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 4
mycols = gg_color_hue(n)

ph1<-ggplot(data = sumdat1, aes(y=MSE, x=Method)) + 
  geom_boxplot(aes(fill=Method))+
  facet_wrap( ~ cens+samp,ncol = 2)+
  labs(y="MSE",x="Method",title = "Scenario 2")+
  theme(legend.position = "none",axis.title.x=element_blank())+
  scale_fill_manual(values=mycols[-4])+
  #theme(axis.text.x=element_text(angle=10,hjust=1,vjust=0.4),plot.title = element_text(hjust = 0.5))
  theme( axis.text = element_text( size = 16 ),
         axis.text.x = element_text( size = 14,angle = 15, hjust=1 ),
         axis.title = element_text( size = 16),
         legend.title = element_text( size = 16 ),
         legend.text=element_text( size = 14 ),
         strip.text = element_text(size = 16),
         plot.title = element_text(size = 18,hjust = 0.5))


###############Chosing the optimal bandwidth
true1<-subset(true1,samp=="n=100")
brier1<-subset(brier1,samp=="n=100")


h.true1<-data.frame(h.grid11,t(aggregate(true1[,(ncol(true1)-hn-2+1):(ncol(true1)-2)],by=list(factor(true1$samp),factor(true1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))
h.brier1<-data.frame(h.grid11,t(aggregate(brier1[,(ncol(brier1)-hn-2+1):(ncol(brier1)-2)],by=list(factor(brier1$samp),factor(brier1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))


colnames(h.true1)<-c("h.grid","n100_c15","n100_c25")
colnames(h.brier1)<-c("h.grid","n100_c15","n100_c25")

h.true1<-h.true1%>% gather(Method, Error, 2:3)
h.brier1<-h.brier1%>% gather(Method, Error, 2:3)


opt.htrue<-data.frame(samp=factor(true1$samp),cens=factor(true1$cens),h.opt=true1$h.opt)
vert.htrue<-data.frame(Method=levels(factor(h.brier1$Method)),h.opt=ddply(opt.htrue,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

opt.hbrier<-data.frame(samp=factor(brier1$samp),cens=factor(brier1$cens),h.opt=brier1$h.opt)
vert.hbrier<-data.frame(Method=levels(factor(h.brier1$Method)),h.opt=ddply(opt.hbrier,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

face_lab<-as_labeller(c("n100_c15"="n=100, CR=15%", "n100_c25"="n=100, CR=25%"))


figure1 <- ggplot(data = h.brier1) + 
  geom_line(aes(x = h.grid , y = Error))+
  geom_vline(data = vert.htrue, aes(xintercept=h.opt,linetype="Brier"),size=0.8)+
  geom_vline(data = vert.hbrier, aes(xintercept=h.opt,linetype="MSE"),size=0.8)+
  #ylim(0.092,0.098)+
  geom_text(x=3, y=30, label="Scatter plot")+
  facet_wrap( ~ Method,ncol=1,labeller = face_lab)+
  labs(y="Brier Score",x="h.grid",title = "Scenario 2")+
  scale_y_continuous(breaks=seq(0.07,0.1,by=0.003))+
  theme(legend.position = "bottom")+
  #geom_text(data = ann_text, mapping= aes(x = 1.5, y=0.0975,label = label),size=5)+
  scale_color_manual(name = "Bandwidth", values = c(MSE = "#F8766D", Brier = "#7CAE00"))+
  scale_linetype_manual(name = "Bandwidth", values = c(MSE = "dashed", Brier = "dotted"))+
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12 ),
         axis.title = element_text( size = 14),
         legend.title = element_text( size = 14 ),
         legend.text=element_text( size = 12 ),
         strip.text = element_text(size = 14),
         plot.title = element_text(size = 14,hjust = 0.5));figure1



############Scenario4

setwd("D:\\UofL\\Project\\Dr.Qi\\Simulations\\Results\\Codes\\CV\\Simu\\standardized_newbr\\Scenario4\\")

true1<-read.table("true1.txt",na.strings = "NA")
brier1<-read.table("brier1.txt",na.strings = "NA")

h.grid3<-seq(0.5,2,by=0.05)

hn<-length(h.grid3)

colnames(true1)<-c("Censoring_rate","h.opt","cond.MSE.mean","KM.MSE.mean","FLCRM.MSE.mean","Cond.MSE.sd","KM.MSE.sd","FLCRM.MSE.sd",paste(rep("h",hn),1:hn,rep(".err",hn),sep = ""))
colnames(brier1)<-colnames(true1)

true1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
true1$cens<-rep(c("15% Censoring","25% Censoring"),200)

brier1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
brier1$cens<-rep(c("15% Censoring","25% Censoring"),200)

sumdat3 <- data.frame(samp=true1$samp,cens=true1$cens,Cond.KM.MSE=true1$cond.MSE.mean,Cond.KM.Brier=brier1$cond.MSE.mean,KM=brier1$KM.MSE.mean,FLCRM=brier1$FLCRM.MSE.mean) %>% gather(Method, MSE, 3:6)
sumdat3$Method<-as.factor(sumdat3$Method)
levels(sumdat3$Method)
sumdat3$Method <- factor(sumdat3$Method , levels=c("Cond.KM.Brier", "Cond.KM.MSE", "KM",  "FLCRM"))


ph3<-ggplot(data = sumdat3, aes(y=MSE, x=Method)) + 
  geom_boxplot(aes(fill=Method))+
  #scale_fill_manual(values=c("#AD0100", "#7cAE00", "#00BFC4")) +
  facet_wrap( ~ cens+samp,ncol = 2)+
  labs(y="MSE",x="Method",title = "Scenario 4")+
  theme(legend.position = "bottom",axis.title.x=element_blank())+
  scale_fill_manual(values=c(mycols))+
  theme( axis.text = element_text( size = 16 ),
         axis.text.x = element_text( size = 14,angle = 15, hjust=1 ),
         axis.title = element_text( size = 16),
         legend.title = element_text( size = 14 ),
         legend.text=element_text( size = 12 ),
         strip.text = element_text(size = 16),
         plot.title = element_text(size = 18,hjust = 0.5));ph3

###############Chosing the optimal bandwidth

true1<-subset(true1,samp=="n=100")
brier1<-subset(brier1,samp=="n=100")

h.true3<-data.frame(h.grid3,t(aggregate(true1[,(ncol(true1)-hn-2+1):(ncol(true1)-2)],by=list(factor(true1$samp),factor(true1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))
h.brier3<-data.frame(h.grid3,t(aggregate(brier1[,(ncol(brier1)-hn-2+1):(ncol(brier1)-2)],by=list(factor(brier1$samp),factor(brier1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))


colnames(h.true3)<-c("h.grid","n100_c15","n100_c25")
colnames(h.brier3)<-c("h.grid","n100_c15","n100_c25")

h.true3<-h.true3%>% gather(Method, Error, 2:3)
h.brier3<-h.brier3%>% gather(Method, Error, 2:3)


opt.htrue<-data.frame(samp=factor(true1$samp),cens=factor(true1$cens),h.opt=true1$h.opt)
vert.htrue<-data.frame(Method=levels(factor(h.brier3$Method)),h.opt=ddply(opt.htrue,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

opt.hbrier<-data.frame(samp=factor(brier1$samp),cens=factor(brier1$cens),h.opt=brier1$h.opt)
vert.hbrier<-data.frame(Method=levels(factor(h.brier3$Method)),h.opt=ddply(opt.hbrier,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

face_lab<-as_labeller(c("n100_c15"="n=100, CR=15%", "n100_c25"="n=100, CR=25%"))


figure3 <- ggplot(data = h.brier3) + 
  geom_line(aes(x = h.grid , y = Error))+
  geom_vline(data = vert.htrue, aes(xintercept=h.opt,linetype="Brier"),size=0.8)+
  geom_vline(data = vert.hbrier, aes(xintercept=h.opt,linetype="MSE"),size=0.8)+
  #ylim(0.092,0.098)+
  geom_text(x=3, y=30, label="Scatter plot")+
  facet_wrap( ~ Method,ncol=1,labeller = face_lab)+
  labs(y="Brier Score",x="h.grid",title = "Scenario 4")+
  theme(legend.position = "bottom")+
  #geom_text(data = ann_text, mapping= aes(x = 1.5, y=0.0975,label = label),size=5)+
  scale_color_manual(name = "Bandwidth", values = c(MSE = "#F8766D", Brier = "#7CAE00"))+
  scale_linetype_manual(name = "Bandwidth", values = c(MSE = "dashed", Brier = "dotted"))+
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12 ),
         axis.title = element_text( size = 14),
         legend.title = element_text( size = 14 ),
         legend.text=element_text( size = 12 ),
         strip.text = element_text(size = 14),
         plot.title = element_text(size = 14,hjust = 0.5));figure3


####Scenario3

setwd("D:\\UofL\\Project\\Dr.Qi\\Simulations\\Results\\Codes\\CV\\Simu\\standardized_newbr\\Scenario3\\")


h.grid2<-seq(1,20,by=1)
hn<-length(h.grid2)

true1<-read.table("true1.txt",na.strings = "NA")
brier1<-read.table("brier1.txt",na.strings = "NA")

colnames(true1)<-c("Censoring_rate","h.opt","cond.MSE.mean","KM.MSE.mean","Cond.MSE.sd","KM.MSE.sd",paste(rep("h",hn),1:hn,rep(".err",hn),sep = ""))
colnames(brier1)<-colnames(true1)

true1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
true1$cens<-rep(c("15% Censoring","25% Censoring"),200)

brier1$samp<-rep(c("n=100","n=100","n=400","n=400"),100)
brier1$cens<-rep(c("15% Censoring","25% Censoring"),200)


sumdat2 <- data.frame(samp=true1$samp,cens=true1$cens,Cond.KM.MSE=true1$cond.MSE.mean,Cond.KM.Brier=brier1$cond.MSE.mean,KM=brier1$KM.MSE.mean) %>% gather(Method, MSE, 3:5)


ph2<-ggplot(data = sumdat2, aes(y=MSE, x=Method)) + 
  geom_boxplot(aes(fill=Method))+
  facet_wrap( ~ cens+samp,ncol = 2)+
  labs(y="MSE",x="Method",title = "Scenario 3")+
  theme(legend.position = "none",axis.title.x=element_blank())+
  scale_fill_manual(values=mycols[-4])+
  #theme(axis.text.x=element_text(angle=10,hjust=1,vjust=0.4),plot.title = element_text(hjust = 0.5))
  theme( axis.text = element_text( size = 10 ),
         axis.text.x = element_text( size = 10,angle = 15, hjust=1 ),
         axis.title = element_text( size = 16),
         legend.title = element_text( size = 16 ),
         legend.text=element_text( size = 14 ),
         strip.text = element_text(size = 16),
         plot.title = element_text(size = 18,hjust = 0.5));ph2


###############Chosing the optimal bandwidth

true1<-subset(true1,samp=="n=100")
brier1<-subset(brier1,samp=="n=100")

h.true2<-data.frame(h.grid2,t(aggregate(true1[,(ncol(true1)-hn-2+1):(ncol(true1)-2)],by=list(factor(true1$samp),factor(true1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))
h.brier2<-data.frame(h.grid2,t(aggregate(brier1[,(ncol(brier1)-hn-2+1):(ncol(brier1)-2)],by=list(factor(brier1$samp),factor(brier1$cens)), mean,na.rm=TRUE)[,-c(1,2)]))


colnames(h.true2)<-c("h.grid","n100_c15","n100_c25")
colnames(h.brier2)<-c("h.grid","n100_c15","n100_c25")

h.true2<-h.true2%>% gather(Method, Error, 2:3)
h.brier2<-h.brier2%>% gather(Method, Error, 2:3)

opt.htrue<-data.frame(samp=factor(true1$samp),cens=factor(true1$cens),h.opt=true1$h.opt)
vert.htrue<-data.frame(Method=levels(factor(h.brier2$Method)),h.opt=ddply(opt.htrue,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

opt.hbrier<-data.frame(samp=factor(brier1$samp),cens=factor(brier1$cens),h.opt=brier1$h.opt)
vert.hbrier<-data.frame(Method=levels(factor(h.brier2$Method)),h.opt=ddply(opt.hbrier,c("samp","cens"),summarise,mean=mean(h.opt))[,3])

face_lab<-as_labeller(c("n100_c15"="n=100, CR=15%", "n100_c25"="n=100, CR=25%"))


figure2 <- ggplot(data = h.brier2) + 
  geom_line(aes(x = h.grid , y = Error))+
  geom_vline(data = vert.htrue, aes(xintercept=h.opt,linetype="Brier"),size=0.8)+
  geom_vline(data = vert.hbrier, aes(xintercept=h.opt,linetype="MSE"),size=0.8)+
  #ylim(0.092,0.098)+
  geom_text(x=3, y=30, label="Scatter plot")+
  facet_wrap( ~ Method,ncol=1,labeller = face_lab)+
  labs(y="Brier Score",x="h.grid",title = "Scenario 3")+
  theme(legend.position = "bottom")+
  #geom_text(data = ann_text, mapping= aes(x = 1.5, y=0.0975,label = label),size=5)+
  scale_color_manual(name = "Bandwidth", values = c(MSE = "#F8766D", Brier = "#7CAE00"))+
  scale_linetype_manual(name = "Bandwidth", values = c(MSE = "dashed", Brier = "dotted"))+
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12 ),
         axis.title = element_text( size = 14),
         legend.title = element_text( size = 14 ),
         legend.text=element_text( size = 12 ),
         strip.text = element_text(size = 14),
         plot.title = element_text(size = 14,hjust = 0.5));figure2

####################################################################

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(ph3)

figure <- grid.arrange(arrangeGrob(ph0+theme(axis.text.x=element_blank(),legend.position = "none",axis.title = element_text( size = 14),
                                             strip.text = element_text(size = 12),
                                             axis.text = element_text( size = 10 ),
                                             strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
                                             plot.title = element_text(size = 14)),
                                   ph1+theme(axis.text.x=element_blank(),
                                             strip.text = element_text(size = 12),
                                             axis.text = element_text( size = 10 ),
                                             axis.title=element_blank(),
                                             strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
                                             plot.title = element_text(size = 14)),
                                   ph2+theme(axis.text.x=element_blank(),
                                             strip.text = element_text(size = 12),
                                             axis.text = element_text( size = 10 ),
                                             axis.title = element_text( size = 14),
                                             strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
                                             plot.title = element_text(size = 14)),
                                   ph3+theme(axis.text.x=element_blank(),legend.position = "none",
                                             legend.title = element_text( size = 14 ),
                                             legend.text=element_text( size = 12 ),
                                             axis.text = element_text( size = 10 ),
                                             axis.title=element_blank(),
                                             strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
                                             strip.text = element_text(size = 12),
                                             plot.title = element_text(size = 14)),nrow=2),
                       mylegend,nrow = 2, heights=c(14,1))

theme(
  axis.title = element_text( size = 18),
  legend.title = element_text( size = 24 ),
  legend.text=element_text( size = 22 ),
  strip.text = element_text(size = 18),
  plot.title = element_text(size = 24,hjust = 0.5))


figure <- grid.arrange(arrangeGrob(figure0,figure1,figure2,figure3+
                                     theme(legend.position = "none"),nrow=2),
                       mylegend,nrow = 2, heights=c(10,1))

figure0<-annotate_figure(figure0,top = text_grob("Scenario 1"))
figure1<-annotate_figure(figure1,top = text_grob("Scenario 2"))
figure2<-annotate_figure(figure2,top = text_grob("Scenario 3"))
figure3<-annotate_figure(figure3,top = text_grob("Scenario 4"))


mylegend<-g_legend(figure3)
figure <- grid.arrange(arrangeGrob(figure0+theme(legend.position = "none",axis.title = element_text( size = 16),
                                                 strip.text = element_text(size = 14),
                                                 plot.title = element_text(size = 16)),
                                   figure1+theme(legend.position = "none",axis.title = element_text( size = 16),
                                                 strip.text = element_text(size = 14),
                                                 plot.title = element_text(size = 16)),
                                   figure2+theme(legend.position = "none",axis.title = element_text( size = 16),
                                                 strip.text = element_text(size = 14),
                                                 plot.title = element_text(size = 16)),
                                   figure3+
                                     theme(legend.position = "none",
                                           axis.title = element_text( size = 16),
                                           strip.text = element_text(size = 14),
                                           plot.title = element_text(size = 16)),nrow=2),mylegend,
                       nrow = 2, heights=c(12,1))

annotate_figure(figure,top.left = text_grob("Scenario 1"))



###########Calculating mean MSE for table


sim.mean<-tapply(sumdat$MSE, list(sumdat$Method,sumdat$cens,sumdat$samp), mean)

sim.sd<-tapply(sumdat$MSE, list(sumdat$Method,sumdat$cens,sumdat$samp), sd)

cbind(as.vector(sim.mean[,1,]),as.vector(sim.sd[,1,]),as.vector(sim.mean[,2,]),as.vector(sim.sd[,2,]))


sim.mean<-tapply(sumdat1$MSE, list(sumdat1$Method,sumdat1$cens,sumdat1$samp), mean)

sim.sd<-tapply(sumdat1$MSE, list(sumdat1$Method,sumdat1$cens,sumdat1$samp), sd)

cbind(as.vector(sim.mean[,1,]),as.vector(sim.sd[,1,]),as.vector(sim.mean[,2,]),as.vector(sim.sd[,2,]))


sim.mean<-tapply(sumdat2$MSE, list(sumdat2$Method,sumdat2$cens,sumdat2$samp), mean)

sim.sd<-tapply(sumdat2$MSE, list(sumdat2$Method,sumdat2$cens,sumdat2$samp), sd)

cbind(as.vector(sim.mean[,1,]),as.vector(sim.sd[,1,]),as.vector(sim.mean[,2,]),as.vector(sim.sd[,2,]))



sim.mean<-tapply(sumdat3$MSE, list(sumdat3$Method,sumdat3$cens,sumdat3$samp), mean)

sim.sd<-tapply(sumdat3$MSE, list(sumdat3$Method,sumdat3$cens,sumdat3$samp), sd)

cbind(as.vector(sim.mean[,1,]),as.vector(sim.sd[,1,]),as.vector(sim.mean[,2,]),as.vector(sim.sd[,2,]))

head(sim.mean)
