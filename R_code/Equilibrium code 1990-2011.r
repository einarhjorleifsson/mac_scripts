####
#
#run with R3.0.3


rm(list=ls())

library(FLCore)


path          <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/"
outPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/R code"
plotPath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/exploration/plots/SRwith EqSim/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Data"

recrPeriod<-1990:2011

StkNam=paste("NEA Mac ",min(recrPeriod),"_",max(recrPeriod),sep="")
require(devtools)
#install_github("msy","wgmg")
require(msy)

load(paste(outPath,"/Mac.RData",sep=""))

stock=window(Mac,start=min(recrPeriod),end=max(recrPeriod))

discards.n(stock)=0
discards(stock)=computeDiscards(stock)
catch(stock)=computeCatch(stock)
landings(stock)=catch(stock)
landings.wt(stock)=catch.wt(stock)
landings.n(stock)=catch.n(stock)
stock(stock)=computeStock(stock)
plot(stock)

FIT <- fitModels(stock,
                 nsamp = 1000,
                 runid = StkNam,
                 method = "Buckland")
SRplot(FIT)
savePlot(filename = paste(outPath,StkNam,"S-R"),type = "png")
graphics.off()



Modset<-FIT$fit
names(Modset)<-c("A","B","sigma","mod")
Modset$A[Modset$mod=="segreg"]<-Modset$A[Modset$mod=="segreg"]*Modset$B[Modset$mod=="segreg"]


BHasave=Modset$A[Modset$mod=="bevholt"];BHbsave= Modset$B[Modset$mod=="bevholt"]
#Transform to the model spec used by FLR   
BHb=BHasave/BHbsave
BHa=1/BHbsave
Modset$A[Modset$mod=="bevholt"]<- BHa
Modset$B[Modset$mod=="bevholt"]<- BHb

Modset$mod<-factor(Modset$mod)



par(mfrow=c(3,2))
for (m in levels(Modset$mod))
{
plot(density(Modset$A[Modset$mod==m]),main=paste(m,"A"))
plot(density(Modset$B[Modset$mod==m]),main=paste(m,"B"))
}
savePlot(filename = paste(outPath,"params distrib"),type = "png")
graphics.off()



# write.csv(Modset,file=paste(inPath,"/SRwithEqSim",min(recrPeriod),"_",max(recrPeriod),".csv",sep=""))





# comparison with the output of the SR bayes analysis
 if(1==2)
 { 
        Modset2<-read.csv(file=paste(inPath,"/SRbayes1990_2011.csv",sep=""))[,-1]
        
   
par(mfrow=c(3,2))
for (m in levels(Modset2$mod))
{
plot(density(Modset2$A[Modset2$mod==m]),main=paste(m,"A"))
plot(density(Modset2$B[Modset2$mod==m]),main=paste(m,"B"))
}     
        
        
        
    Modset3<-read.csv(file=paste(inPath,"/SRIndivFit.csv",sep=""))[,-1]
        
   
par(mfrow=c(3,2))
for (m in levels(Modset3$mod))
{
plot(density(Modset3$A[Modset3$mod==m]),main=paste(m,"A"))
plot(density(Modset3$B[Modset3$mod==m]),main=paste(m,"B"))
}    
    
   dim(Modset3)
    
    
    
    
    
    
        
        mods<-Modset
        mods$meth<-"EqSim"
        mods2<-Modset2
        mods2$meth<-"SR Bayes"
        mods2<-mods2[,c(2,3,4,1,5)]
        modss<-rbind(mods,mods2)
        
        histogram(A~1|meth,data=modss[modss$mod=="segreg",])
        densityplot(A~1|meth,data=modss[modss$mod=="segreg",])
        
     
        HS<-Modset[Modset$mod=="segreg",]
        HS$meth<-"EqSim"
        HS2<-Modset2[Modset2$mod=="segreg",]
        HS2$meth<-"SR Bayes"
        HS2<-HS2[,c(2,3,4,1,5)]
        HSs<-rbind(HS,HS2)
     
     
     
     
        summary(HS);summary(HS2)
        
        
        
        BH<-Modset[Modset$mod=="bevholt",]
        BH2<-Modset2[Modset2$mod=="bevholt",]
        summary(BH);summary(BH2)
        
        
        RK<-Modset[Modset$mod=="ricker",]
        RK2<-Modset2[Modset2$mod=="ricker",]
        summary(RK);summary(RK2)
  }
















