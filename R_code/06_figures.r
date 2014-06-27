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
library(FLAssess)
library(FLSAM)
library(MASS)
library(msm)

wine <- F

path          <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/R code/"
outPath       <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/"
plotPath       <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/Plots input/"
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



#- Load objects
load(file=paste(outPath,"Mac.RData",            sep=""))
load(file=paste(outPath,"Macctrl.RData",        sep=""))

load(file=paste(outPath,"biol.RData",           sep=""))
load(file=paste(outPath,"fishery.RData",        sep=""))
  load(file=paste(outPath,"propN.RData",        sep=""))
  load(file=paste(outPath,"propWt.RData",       sep=""))
  load(file=paste(outPath,"ctch.RData",         sep=""))
load(file=paste(outPath,"surveys.RData",        sep=""))
  load(file=paste(outPath,"surv.RData",         sep=""))
  load(file=paste(outPath,"surveyQ.RData",      sep=""))
  load(file=paste(outPath,"surveyK.RData",      sep=""))
load(file=paste(outPath,"stocks.RData",         sep=""))
load(file=paste(outPath,"settings.RData",       sep=""))
load(file=paste(outPath,"resNFinal.RData",      sep=""))
load(file=paste(outPath,"SRmod.RData",      sep=""))
load(file=paste(outPath,"resFFinal.RData",      sep=""))
for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])


scen          <- c("LTMP")              # 
opt           <- 1                      # for multiple scenario combinations a counter
TACvarlim     <- T                      # whether or not to apply the 20% limit on TAC change
Fvarlim       <- T                      # whether or not to apply the 10% limit on Fbar change
BBscen        <- "AlternateBank"        # banking borrowing options :
                                                     # "Banking"          : always bank 
                                                     # "Borrowing"        : always borrow
                                                     # "AlternateBank"    : bank first and then alternate
                                                     # "AlternateBorrow"  : borrow first and then alternate
                                                     # "MinVar"           : use BB to minise TAC variability
                                                     
LastRecOp     <- "RCT3"                  # option to replace the last estimated recruitment : "SAM", "geom", "RCT3"
                                                     # "SAM"  = don't overwrite the SAM estimmate
                                                     # "SAM"  = don't overwrite the SAM estimmate
                                                     # "geom" = replace by geomean 1990/(TaY-1)
                                                     # "RCT3" = replace by RCT3 output
mpOptions<-list(opt=opt,TACvarlim=TACvarlim,Fvarlim=Fvarlim,BBscen=BBscen)

sc.name<-paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep="")

outPath<-paste(outPath,sc.name,"/",sep="")

scen<-paste(scen,opt,sep="")

load(file=paste(outPath,"/",scen,"_Finalf.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_Finalbiol.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_Finalstocks.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_Finalpercievedstocks.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_FinalTAC.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_FinalfSTF.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_Finalfishery.RData",           sep=""))
load(file=paste(outPath,"/",scen,"_FinalmpPoints.RData",           sep=""))
source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))

RecType<-settings$RecType

##-------------------------------------------------------------------------------
## Figures on stochastic variables
##-------------------------------------------------------------------------------

ca    <- 1.1
cl    <- 1.3
fonts <- 2

  #-------------------------------------------------------------------------------
  # 0): Figure of biol future weights and recruitment
  #-------------------------------------------------------------------------------
plotPath<-paste(outPath,"Plots input/",sep="")
dir.create(plotPath)


    #- Recruitment lognormal default fit
hist(biol@n[1,projPeriod],breaks=100,xlab="Recruitment",ylab="Simulated frequency",main="",cex.lab=cl,cex.axis=ca,font=fonts)
abline(v=exp(mean(log(biol@n[1,recrPeriod]))),col="red",lty=3,lwd=3)
legend("topright",legend=c("Simulated recruitment","Geometric mean recruitment"),lty=c(0,2),pch=c(22,0),
       pt.bg="white",pt.cex=c(1,0),box.lty=0,col=c("black","red"),lwd=c(0,3))
savePlot(file=paste(plotPath,"input_recruitment.png",sep=""),type="png")

    #- Recruitment cummulative distribution
plot(x=sort(c(rec(Mac[,ac(recrPeriod)]))),y=seq(0,1,length.out=length(c(rec(Mac[,ac(recrPeriod)])))+2)[-c(1,length(c(rec(Mac[,ac(recrPeriod)])))+2)],
     ylim=c(0,1),col=2,pch=15,xlab="Recruitment",ylab="Probability",las=1,xlim=range(c(biol[1,projPeriod]@n)),yaxs="i")
points(x=sort(c(biol@n[1,ac(projPeriod)])),y=seq(0,1,length.out=length(sort(c(biol@n[1,ac(projPeriod)])))+2)[c(-1,-length(sort(c(biol@n[1,ac(projPeriod)]))))],pch=19,cex=0.3)
legend("topleft",legend=c("Simulated recruitment","Observed recruitment"),lty=c(0,0),pch=c(16,15),
       pt.cex=c(0.5,1),box.lty=0,col=c("black","red"),lwd=c(0,0))
savePlot(file=paste(plotPath,"input_recruitmentDistri.png",sep=""),type="png")

    #- SR plot with simulated data
Recs<-c(biol@n[1,ac(projPeriod)]@.Data)
SSBs<-c(ssbb(biol,f,stocks)[,ac(projPeriod)]@.Data)/1e6
RecO<-c(rec(Mac[,ac(recrPeriod)]))
SSBO<-c(ssb(Mac[,ac(recrPeriod)]))/1e6


png(file=paste(plotPath,"input_simulated recruitments.png",sep=""),width=8,height=6,units="in",res=300)
plot(SSBs,Recs,xlim=c(0,8),xlab="SSB (10^6t)",ylab="Rec (thousands)",cex=0.4,pch=19,main=paste("RecType=",RecType,"Ftarget=",mpPoints$Ftarget))
seqMin <- 0.5
seqStep <- 0.5
Xmax <-8
up=seq(seqMin,Xmax,seqStep)
lw=seq(seqMin,Xmax,seqStep)
md=seq(seqMin,Xmax,seqStep)
ssb=seq(seqMin,Xmax,seqStep)
loopNum <- length(ssb)    
for (j in 1:loopNum){     
  up[j]=quantile(Recs[((SSBs>up[j])*(SSBs<up[j+1]))>0],probs=.95,na.rm=TRUE)
  lw[j]=quantile(Recs[((SSBs>lw[j])&(SSBs<lw[j+1]))>0],probs=.05,na.rm=TRUE)
  md[j]=quantile(Recs[((SSBs>md[j])&(SSBs<md[j+1]))>0],probs=.5,na.rm=TRUE)
  }
lines(ssb,md[1:loopNum],col="green",lwd=3)
lines(ssb,up[1:loopNum],col=4,lwd=3)
lines(ssb,lw[1:loopNum],col=4,lwd=3)
points(SSBO,RecO,col="red",cex=1,pch=19)
dev.off()


#
#
#RecO<-c(rec(Mac[,ac(recrPeriod)]))
#SSBO<-c(ssb(Mac[,ac(recrPeriod)]))/1e6
#Rec<-c(rec(Mac))
#SSB<-c(ssb(Mac))/1e6
#
#plot(SSB,Rec,xlim=c(0,6),ylim=c(0,1.2e7),xlab="SSB (10^6t)",ylab="Rec (thousands)",cex=0.4,pch=19,main="")
#points(SSBO,RecO,col="red",cex=1,pch=19)
#





    #- Numbers-at-age
xyplot(data ~ age | as.factor(year),as.data.frame(biol@n),groups=iter,type="p",col=1,pch=19,cex=0.3)
savePlot(file=paste(plotPath,"input_natageXY.png",sep=""),type="png"); dev.off()
    
par(mfrow=c(7,8),mar=c(0.3,0.1,0.3,0.1),oma=c(6,6,2,2),las=1)
  for(iYr in 1980:2012){
    if(iYr %in% seq(1980,2007,8)) boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),xaxt="n",xlab="",ylab="",ylim=c(0,rev(pretty(max(as.data.frame(biol@n)$data,na.rm=T)))[1]))
    if(iYr %in% 2009:2010)        boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),yaxt="n",xlab="",ylab="",ylim=c(0,rev(pretty(max(as.data.frame(biol@n)$data,na.rm=T)))[1]))
    if(iYr == 2008)               boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),xlab="",ylab="",ylim=c(0,rev(pretty(max(as.data.frame(biol@n)$data,na.rm=T)))[1]))
    if(!iYr %in% c(seq(1980,2007,8),2008:2010)) boxplot(data ~ as.factor(age),data=as.data.frame(biol@n[,ac(iYr)]),xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,rev(pretty(max(as.data.frame(biol@n)$data,na.rm=T)))[1]))
    text(x=9,y=rev(pretty(max(as.data.frame(biol@n)$data,na.rm=T)))[1],pos=1,labels=iYr,font=2,cex=1)
}
savePlot(file=paste(plotPath,"input_natage.png",sep=""),type="png"); dev.off()


  #- Weight at age
boxplot(as.data.frame(biol@wt[,projPeriod])$data~as.factor(as.data.frame(biol@wt[,projPeriod])$age),
        boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
        xlab="Age",ylab="Weight in the stock(kg)",cex.axis=ca,cex.lab=cl,font=fonts)
points(as.data.frame(iter(Mac@stock.wt[,ac(2001:2012)],1))$data~as.factor(as.data.frame(iter(Mac@stock.wt[,ac(2001:2012)],1))$age),pch=19,cex=0.5,col="red")
legend("bottomright",legend=c("Simulated weights","Observed weights"),col=c("black","red"),
       pch=c(22,19),
       pt.bg="grey90",pt.cex=c(1,1),box.lty=0)
savePlot(file=paste(plotPath,"input_wtatage.png",sep=""),type="png"); dev.off()

xyplot(data ~ year | as.factor(age),data=subset(subset(as.data.frame(biol@wt[ac(1:6)]),year %in% 2001:futureMaxYr),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_wtscenarios.png",sep=""),type="png"); dev.off()

  #- Maturity at age
boxplot(as.data.frame(biol@fec[,projPeriod])$data~as.factor(as.data.frame(biol@fec[,projPeriod])$age),
        boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
        xlab="Age",ylab="Proportion mature",cex.axis=ca,cex.lab=cl,font=fonts)
points(as.data.frame(iter(Mac@mat[,ac(2001:2010)],1))$data~as.factor(as.data.frame(iter(Mac@mat[,ac(2001:2010)],1))$age),pch=19,cex=0.5,col="red")
legend("bottomright",legend=c("Simulated maturity","Observed maturity"),col=c("black","red"),
       pch=c(22,19),
       pt.bg="grey90",pt.cex=c(1,1),box.lty=0)
savePlot(file=paste(plotPath,"input_matatage.png",sep=""),type="png"); dev.off()

xyplot(data ~ year | as.factor(age),data=subset(subset(as.data.frame(biol@fec),year %in% 2001:2022),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(list(...)$y)))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_matscenarios.png",sep=""),type="png"); dev.off()



  #-------------------------------------------------------------------------------
  # 1): Figure of fisheries future landing weights and selectivity
  #-------------------------------------------------------------------------------

  #- Landings weight

for(iFsh in dimnames(fishery@landings.wt)$unit){
  boxplot(as.data.frame(fishery@landings.wt[,projPeriod,iFsh])$data~as.factor(as.data.frame(fishery@landings.wt[,projPeriod,iFsh])$age),
          boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
          xlab="Age",ylab="Weight (kg)",cex.lab=cl,cex.axis=ca,font=fonts,main=paste("Fleet",iFsh))
  points(as.data.frame(iter(fishery@landings.wt[,ac(2001:2010),iFsh],1))$data~as.factor(as.data.frame(iter(fishery@landings.wt[,ac(2001:2010),iFsh],1))$age),pch=19,cex=0.5,col="red")
  legend("bottomright",legend=c("Simulated weights","Observed weights"),col=c("black","red"),
         pch=c(22,19),
         pt.bg="grey90",pt.cex=c(1,1),box.lty=0)
}
savePlot(file=paste(plotPath,"input_fish_wtatage.png",sep=""),type="png"); dev.off()

    #- Selectivity

  boxplot(as.data.frame(landings.sel(fishery)[,projPeriod])$data~as.factor(as.data.frame(landings.sel(fishery)[,projPeriod])$age),
          boxwex=0.5,col="grey90",medlty=0,medlwd=1,medpch=19,medcol="black",boxlty=1,outlty=0,outcex=0.5,outpch=1,staplelty=1,
          xlab="Age",ylab="Selectivity",cex.lab=cl,cex.axis=ca,font=fonts,main=paste("Fleet",iFsh))
  lines(as.data.frame(apply(landings.sel(fishery)[,projPeriod],1,median))$data~as.factor(as.data.frame(apply(landings.sel(fishery)[,projPeriod],1,median))$age),
        lwd=2,lty=2)

  lines(as.data.frame(iter(landings.sel(fishery)[,ac(histMaxYr)],1))$data~as.factor(as.data.frame(iter(landings.sel(fishery)[,ac(histMaxYr)],1))$age),
        col="red",lty=1,lwd=2)

savePlot(file=paste(plotPath,"input_fish_selatage.png",sep=""),type="png"); dev.off()

xyplot(data ~ age | as.factor(year),data=subset(subset(as.data.frame(landings.sel(fishery)[,,1]),year %in% 2001:futureMaxYr),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_selAscenarios.png",sep=""),type="png"); dev.off()

  #-------------------------------------------------------------------------------
  # 2): Figure of survey residual pattern
  #-------------------------------------------------------------------------------


xyplot(data ~ year,subset(subset(as.data.frame(Rindex[1,]),year %in% 1998:2012),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_Rindex_errors.png",sep=""),type="png"); dev.off()



  #-------------------------------------------------------------------------------
  # 3): Figures of the stock
  #-------------------------------------------------------------------------------

plot(trim(stocks,year=an(histPeriod)))
savePlot(file=paste(plotPath,"input_stocks.png",sep=""),type="png"); dev.off()


devSSB<-log(ssbb(biol,f,stocks)/ssb(stockstore))
devSSB<-devSSB[,1:(dim(devSSB)[2]-1)]

rhoSSBdev<-rep(NA,nits)
for (i in 1:nits) rhoSSBdev[i]<- an(unlist(acf(iter(devSSB,i),plot=F)[1])[1])


xyplot(data ~ year ,data=subset(subset(as.data.frame(devSSB),year %in% 2012:2032),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_SSB_errorscenario.png",sep=""),type="png"); dev.off()





devFbar<-log(quantMeans(f[ac(4:8),])/quantMeans(harvest(stockstore)[ac(4:8),]))
devFbar<-devFbar[,1:(dim(devFbar)[2]-1)]
rhoFdev<-rep(NA,nits)
for (i in 1:nits) rhoFdev[i]<- an(unlist(acf(iter(devFbar,i),plot=F)[1])[1])


xyplot(data ~ year ,data=subset(subset(as.data.frame(devFbar),year %in% 2012:2032),iter %in% 1:5),type="l",group=iter,
       prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
       scales=list(alternating=1,y=list(relation="free",rot=0)))
savePlot(file=paste(plotPath,"input_Fbar_errorscenario.png",sep=""),type="png"); dev.off()

rhos<-data.frame(quel=c(rep("rhoSSBdev",nits),rep("rhoFdev",nits)),rho=c(rhoSSBdev,rhoFdev))
histogram(1~rho|quel, data=rhos,col="grey")
savePlot(file=paste(plotPath,"F_and_SSB_dev_rhos.png",sep=""),type="png"); dev.off()


#- Assessment stock.n error
xyplot(data ~ year|as.factor(age),data=as.data.frame(devN),pch=19,col="black",cex=0.4,xlab="Years",ylab="Error (log ratio)",main="Assessment numbers-at-age error")
savePlot(file=paste(plotPath,"input_stock_errornatage.png",sep=""),type="png"); dev.off()

#- Assessment harvest error
xyplot(data ~ year|as.factor(age),data=as.data.frame(devF),pch=19,col="black",cex=0.4,xlab="Years",ylab="Error (log ratio)",main="Assessment f-at-age error")
savePlot(file=paste(plotPath,"input_stock_errorfatage.png",sep=""),type="png"); dev.off()

#- Catch number residuals
xyplot(data ~ year|as.factor(age),data=as.data.frame(catchResids),pch=19,col="black",cex=0.4,xlab="Years",ylab="Error multiplier",main="Catch numbers-at-age error")
savePlot(file=paste(plotPath,"input_stock_errorcatchatage.png",sep=""),type="png"); dev.off()

#-------------------------------------------------------------------------------
# Figures on results
#-------------------------------------------------------------------------------

plotPath2       <- plotPath<-paste(outPath,"Plot results/",sep="")
dir.create(plotPath2)

rSSBp <- apply(ssbb(biol,unitSums(f),stocks)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rSSBs <- apply(ssb(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFp   <- apply(apply(unitSums(f)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFs   <- apply(fbar(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFstf <- apply(apply(unitSums(fSTF)[ac(4:8),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rLandf<- apply(computeLandings(fishery)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rLands<- apply(computeLandings(stockstore)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rTAC  <- apply(TAC@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)

  #- Plot settings
yrs   <- 2001:2031
cl    <- 1.2
ca    <- 1.1
fonts <- 2


#-------------------------------------------------------------------------------
# 3): Plot results of SSB stock and SSB pop
#-------------------------------------------------------------------------------
par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rSSBp[,,ac(yrs),,,],rSSBs[,,ac(yrs),,,])),na.rm=T)[2])
  #---------
  #- Biology
  #---------  
plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
mtext(text="SSB (thousand tonnes)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Blim & Bpa
abline(h=1.84e6,col="blue",lwd=2,lty=2);
abline(h=2.36e6,col="darkgreen",lwd=2,lty=2);
mtext(text="Blim",side=4,at=1.84e6,las=1,cex=0.8,col="blue",line=0.5,font=fonts)
mtext(text="Bpa",side=4,at=2.36e6,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rSSBp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rSSBs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rSSBs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rSSBs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("SSB assessed stock","SSB true population"),
       col=c("red","black"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,scen,"_SSB.png",sep=""),type="png");dev.off()
#-------------------------------------------------------------------------------
# 4): Plot results of F stock and true F
#-------------------------------------------------------------------------------
par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rFp[,,ac(yrs),,,],rFs[,,ac(yrs),,,],rFstf[,,ac(yrs),,,])),na.rm=T)[2])

  #---------
  #- Biology
  #---------
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text="Fishing mortality (ages 4-8)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Fpa
abline(h=0.26,col="darkgreen",lwd=2,lty=2);
mtext(text="Fpa",side=4,at=0.26,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rFs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rFs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rFs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- STF
  #---------
lines(rFstf["50%",,ac(2014:rev(projPeriod)[2]),,,]~ac(2014:rev(projPeriod)[2]),lty=1,lwd=2,col="blue")
lines(rFstf["5%",,ac(2014:rev(projPeriod)[2]),,,]~ac(2014:rev(projPeriod)[2]),lty=3,col="blue")
lines(rFstf["95%",,ac(2014:rev(projPeriod)[2]),,,]~ac(2014:rev(projPeriod)[2]),lty=3,col="blue")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("F assessed stock","F true population","F short term forecast"),
       col=c("red","black","blue"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,scen,"_F.png",sep=""),type="png");dev.off()
#-------------------------------------------------------------------------------
# 5): Plot results of landings by fleet & TAC on top
#-------------------------------------------------------------------------------


  #---------
  #- Fishery
  #---------
  par(mar=c(5.1,5.1,4.1,2.1))
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),iFsh,,],rTAC[,,ac(yrs),iFsh,,])),na.rm=T)[2])

    #- Landings
  plot(rLandf["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="Year",ylab="Landings (thousand tonnes)",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
  grid(); box()
  lines(rLandf["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rLandf["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rLandf["95%",,ac(yrs),iFsh,,]~yrs,lty=3)
    #- TAC
  lines(rTAC["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2,col="red")
  lines(rTAC["5%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")
  lines(rTAC["95%",,ac(yrs),iFsh,,]~yrs,lty=3,col="red")



legend("bottomright",legend=c("Landings","Advice TAC"),
       col=c("black","red"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,scen,"_CatchTAC.png",sep=""),type="png");dev.off()

#-------------------------------------------------------------------------------
# 6): Plot results total landings fishery and stock
#-------------------------------------------------------------------------------
par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),,,],rLands[,,ac(yrs),,,])),na.rm=T)[2])

  #---------
  #- Biology
  #---------
plot((rLandf["50%",,ac(yrs),,,])~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text="Landings (tonnes)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
lines((rLandf["50%",,ac(yrs),,,])~yrs,lty=1,lwd=2)
lines((rLandf["5%",,ac(yrs),,,])~yrs,lty=3)
lines((rLandf["95%",,ac(yrs),,,])~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rLands["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rLands["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rLands["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("Catch assessed stock","Catch true population"),
       col=c("red","black"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(plotPath2,scen,"_Catch.png",sep=""),type="png");dev.off()


#-------------------------------------------------------------------------------
# 7): Plot results of TAC
#-------------------------------------------------------------------------------

par(mar=c(2.1,2.1,2.1,2.1),oma=c(3,3,2,0))

  #---------
  #- Fishery
  #---------
iFsh<-1
  xrange  <- range(yrs)
  yrange  <- c(0,range(pretty(rTAC[,,ac(yrs),iFsh,,]),na.rm=T)[2])

  plot(rTAC["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
  mtext(text="TAC (thousand tonnes)",side=2,at=0.5,outer=T,cex=cl,line=1.5)
  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
  grid(); box()
  lines(rTAC["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
  lines(rTAC["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
  lines(rTAC["95%",,ac(yrs),iFsh,,]~yrs,lty=3)


savePlot(file=paste(plotPath2,scen,"_TAC.png",sep=""),type="png");dev.off()


#-------------------------------------------------------------------------------
# 8): Plot trajectories of ssb(biol), rec(biol), f(biol), TAC(A)
#-------------------------------------------------------------------------------

par(mfrow=c(2,2),mar=c(3,3,3,3),oma=c(3.1,3.1,1,1))

#- Plot ssb based on biol
xrange  <- range(yrs)
yrange  <- c(0,pretty(range(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)]),stocks[,ac(yrs)]),1:10))/1000),1)[3])
plot(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)]),stocks[,ac(yrs)]),1))/1000~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="Spawning Stock Biomass (kt)",outer=F,line=4,font=fonts,cex=ca)
grid(); box()
for(i in 1:10) lines(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)]),stocks[,ac(yrs)]),i))/1000~yrs,col=i,lwd=2)

#- Plot rec based on biol
xrange  <- range(yrs)
yrange  <- c(0,rev(pretty(range(c(iter(biol@n[1,ac(yrs)],1:10))/1e6),1))[1])
plot(c(iter(biol@n[1,ac(yrs)],1))/1e6~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="Recruitment (millions)",outer=F,line=4,font=fonts,cex=ca)
grid();box()
for(i in 1:10) lines(c(iter(biol@n[1,ac(yrs)],i))/1e6~yrs,col=i,lwd=2)

#- f-ages 2-6
xrange  <- range(yrs)
yrange  <- c(0,pretty(range(c(iter(quantMeans(unitSums(f[ac(4:8),ac(yrs)])),1:10))),1)[3])
plot(c(iter(quantMeans(unitSums(f[ac(4:8),ac(yrs)])),1))~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="Fishing mortality",outer=F,line=4,font=fonts,cex=ca)
grid();box()
for(i in 1:10) lines(c(iter(quantMeans(unitSums(f[ac(4:8),ac(yrs)])),i))~yrs,col=i,lwd=2)

#- TAC of fleet A
xrange  <- range(yrs)
yrange  <- c(0,pretty(range(c(iter(TAC[,ac(yrs),"A"]/1000,1:10))),1)[3])
plot(c(iter(TAC[,ac(yrs),"A"]/1000,1))~yrs,type="l",xlab="Years",ylab="",
     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
mtext(side=2,at=yrange[2]/2,text="TAC (kt)",outer=F,line=4,font=fonts,cex=ca)
grid();box()
for(i in 1:10) lines(c(iter(TAC[,ac(yrs),"A"]/1000,i))~yrs,col=i,lwd=2)

savePlot(file=paste(plotPath2,scen,"_iterations.png",sep=""),type="png");dev.off()




#-------------------------------------------------------------------------------
# 9): Report figures
#-------------------------------------------------------------------------------

yrs   <- 2005:an(rev(projPeriod)[2])
cl    <- 1.1
ca    <- 1
fonts <- 1
yrangeSSB <- c(0,6e6)
yrangeLan <- c(0,2e6)
#yrangeLan2<- c(0,2.5e4)
yrangeF   <- c(0,0.5)


par(mfrow=c(3,1),oma=c(6,6,2,3),mar=c(1,0,0,0))

  #---------
  #- Landings
  #---------

xrange  <- range(yrs)
yrange  <- yrangeLan

plot(rLandf["50%",,ac(yrs),1,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
mtext(text=expression(paste("Landings (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rLandf["50%",,ac(yrs),1,,]~yrs,lty=1,lwd=2)
lines(rLandf["5%",,ac(yrs),1,,]~yrs,lty=3,lwd=1)
lines(rLandf["95%",,ac(yrs),1,,]~yrs,lty=3,lwd=1)
par(new=T)
yrange  <- yrangeLan2
plot(rLandf["50%",,ac(yrs),2,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n",col="red")
axis(4,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts,col="red",col.ticks="red",col.axis="red")
lines(rLandf["50%",,ac(yrs),2,,]~yrs,lty=1,lwd=2,col="red")
lines(rLandf["5%",,ac(yrs),2,,]~yrs,lty=3,lwd=1,col="red")
lines(rLandf["95%",,ac(yrs),2,,]~yrs,lty=3,lwd=1,col="red")
legend("bottomleft",legend=c("Fleet A","Fleet B"),col=c("black","red"),lwd=2,lty=1,box.lty=0)
text(x=xrange[1],y=yrange[2],pos=1,labels="(A)",font=fonts,cex=cl)

  #------------------
  # True F
  #------------------

xrange  <- range(yrs)
yrange  <- yrangeF
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text=expression(paste(F[2-6]," (",year^-1,")",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)

  #-Reference level Fpa
abline(h=0.26,col="darkgreen",lwd=1,lty=2);
mtext(text="Fpa",side=4,at=0.26,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
text(x=xrange[1],y=yrange[2],pos=1,labels="(B)",font=fonts,cex=cl)

  #------------------
  # True SSB
  #------------------

xrange  <- range(yrs)
yrange  <- yrangeSSB
plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),cex=ca,font=fonts)
mtext(text=expression(paste("SSB (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rSSBp["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
  #-Reference level Blim & Bpa
abline(h=1.84e6,col="blue",lwd=1,lty=2);
abline(h=2.36e6,col="darkgreen",lwd=1,lty=2);
mtext(text="Blim",side=4,at=1.84e6,las=1,cex=0.65,col="blue",line=0.5,font=fonts)
mtext(text="Bpa",side=4,at=2.36e6,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
text(x=xrange[1],y=yrange[2],pos=1,labels="(C)",font=fonts,cex=cl)



  #- Labels x-axis
mtext(text=expression(Years),side=1,at=(xrange[2]-xrange[1])/2+xrange[1],outer=F,cex=cl,line=4,font=fonts)
savePlot(paste(plotPath2,scen,"Truth.png",sep=""),type="png")
