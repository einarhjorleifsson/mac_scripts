#-------------------------------------------------------------------------------
# WKHERMPII
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 03-Oct-2011
#
# Build for R2.8.1
#-------------------------------------------------------------------------------

rm(list=ls())
library(FLCore)
#library(FLAssess)
#library(FLSAM)
library(MASS)
library(msm)

wine <- F

path          <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/R code/"
outPath       <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/"


if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/","/media/n/",path)
  inPath      <- sub("W:/","/media/n/",inPath)
  codePath    <- sub("W:/","/media/n/",codePath)
  outPath     <- sub("W:/","/media/n/",outPath)
}


home<-F
if(home)
{
path          <- "D://MSE/"
inPath        <- "D://MSE/Data/"
codePath      <- "D://MSE/R code/"
outPath       <- "D://MSE/Results/"
}
             

source(paste(codePath,"functions.r",            sep=""))
load(file=paste(outPath,"stocks.RData",         sep=""))
 outPath2      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/equilibrium/true/"
 scen          <- "equilTrue"              # 
 opt           <- ""                      # for multiple scenario combinations a counter


SSBs<-c()
Recs<-c()
Yields<-c()
trueF<-c()
fs<- seq(0.0,0.6,0.01)

for (fequ in fs)
{              
              mpPoints      <- list(Fequ=fequ)
cat("! these results are projecting only 20 years head !","\n")      
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalbiol.RData",        sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalfishery.RData",     sep=""))
#load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalpercievedstocks.RData",      sep=""))
load(file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalf.RData",           sep=""))

tp<-ssbb(biol,f,stocks)
tp<-c(yearMeans(tp[,(dim(tp)[2]-5):dim(tp)[2]])@.Data)
SSBs<-c(SSBs,tp)

tp<-n(biol)[ac(1),]
tp<-c(yearMeans(tp[,(dim(tp)[2]-5):dim(tp)[2]])@.Data)
Recs<-c(Recs,tp)

tp<-quantSums(fishery@landings.n*fishery@landings.wt)
tp<-c(yearMeans(tp[,(dim(tp)[2]-5):dim(tp)[2]])@.Data)
Yields<-c(Yields,tp)

tp<-quantMeans(f[ac(4:8),])
tp<-c(yearMeans(tp[,(dim(tp)[2]-5):dim(tp)[2]])@.Data)
trueF<-c(trueF,tp)


}


nits<-dim(f)[6]
iter<-rep(1:nits,length(fs))

res<-data.frame(trueF=trueF,SSB=SSBs,Yields=Yields,Rec=Recs,iter=iter)


# plot as done by Einar
png(file=paste(outPath2,"equilibrium.png",sep=""),width=8,height=8,units="in",res=200)
par(mfrow=c(2,2))
# plot for the recs
q95<-aggregate(Rec~trueF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(Rec~trueF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(Rec~trueF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(Rec~trueF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(Rec~trueF,data=res,function(x) {quantile(x,0.05)})

plot(q50$trueF,q50$Rec,type="l",col="red", lwd=2, xlab="true Fbar",ylab="Rec",main="equilibrium true F imposed",ylim=range(c(q05,q95)))
polygon(c(q25$trueF,rev(q75$trueF)),c(q25$Rec,rev(q75$Rec)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$trueF,rev(q95$trueF)),c(q05$Rec,rev(q95$Rec)),col=rgb(1,0,0,0.3),border=F)
lines(q50$trueF,q50$Rec,lwd=2)

# plot for the ssb
q95<-aggregate(SSB~trueF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(SSB~trueF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(SSB~trueF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(SSB~trueF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(SSB~trueF,data=res,function(x) {quantile(x,0.05)})

plot(q50$trueF,q50$SSB,type="l",col="red", lwd=2, xlab="true Fbar",ylab="SSB",main="equilibrium true F imposed",ylim=range(c(q05,q95)))
polygon(c(q25$trueF,rev(q75$trueF)),c(q25$SSB,rev(q75$SSB)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$trueF,rev(q95$trueF)),c(q05$SSB,rev(q95$SSB)),col=rgb(1,0,0,0.3),border=F)
lines(q50$trueF,q50$SSB,lwd=2)

# plot for the true F
q95<-aggregate(trueF~trueF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(trueF~trueF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(trueF~trueF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(trueF~trueF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(trueF~trueF,data=res,function(x) {quantile(x,0.05)})

plot(q50$trueF,q50$trueF,type="l",col="red", lwd=2, xlab="true Fbar",ylab="trueF",main="equilibrium true F imposed",ylim=range(c(q05,q95)))
polygon(c(q25$trueF,rev(q75$trueF)),c(q25$trueF,rev(q75$trueF)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$trueF,rev(q95$trueF)),c(q05$trueF,rev(q95$trueF)),col=rgb(1,0,0,0.3),border=F)
lines(q50$trueF,q50$trueF,lwd=2)


# plot for the yields
q95<-aggregate(Yields~trueF,data=res,function(x) {quantile(x,0.95)})
q75<-aggregate(Yields~trueF,data=res,function(x) {quantile(x,0.75)})
q50<-aggregate(Yields~trueF,data=res,function(x) {quantile(x,0.50)})
q25<-aggregate(Yields~trueF,data=res,function(x) {quantile(x,0.25)})
q05<-aggregate(Yields~trueF,data=res,function(x) {quantile(x,0.05)})

plot(q50$trueF,q50$Yields,type="l",col="red", lwd=2, xlab="true Fbar",ylab="Yields",main="equilibrium true F imposed",ylim=range(c(q05,q95)))
polygon(c(q25$trueF,rev(q75$trueF)),c(q25$Yields,rev(q75$Yields)),col=rgb(1,0,0,0.6),border=F)
polygon(c(q05$trueF,rev(q95$trueF)),c(q05$Yields,rev(q95$Yields)),col=rgb(1,0,0,0.3),border=F)
lines(q50$trueF,q50$Yields,lwd=2)

dev.off()


MSY<-data.frame(Bmsy=NA,Fmsy=NA,Msy=NA,iter=NA)
for (its in 1:nits)
{
sub<-res[res$iter==its,]
msy<-rev(order(sub$Yields))[1]
MSY<-rbind(MSY,data.frame(Bmsy=sub$SSB[msy],Fmsy=sub$trueF[msy],Msy=sub$Yields[msy],iter=its))
}
MSY<-MSY[-1,]


png(file=paste(outPath2,"iterations msy estimates.png",sep=""),width=12,height=6,units="in",res=200)
par(mfrow=c(1,3))
hist(MSY$Bmsy,nclass=12,density=15,main="Bmsy",xlab="",xlim=c(0,10e6))
abline(v=mean(MSY$Bmsy),lty="dashed",lwd=2)
abline(v=median(MSY$Bmsy),lwd=2)
hist(MSY$Fmsy,nclass=12,density=15,main="Fmsy",xlab="")
abline(v=mean(MSY$Fmsy),lty="dashed",lwd=2)
abline(v=median(MSY$Fmsy),lwd=2)

hist(MSY$Msy,nclass=12,density=15,main="MSY",xlab="")
abline(v=mean(MSY$Msy),lty="dashed",lwd=2)
abline(v=median(MSY$Msy),lwd=2)

dev.off()