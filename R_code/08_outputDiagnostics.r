#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 02-Sep-2012
#
# Build for R2.13.2
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(MASS)
library(msm)
wine <- F

library(FLCore)
library(PBSadmb)
library(lattice)
library(MASS)

ac<-function(x) {return(as.character(x))}
an<-function(x) {return(as.numeric(x))}


path          <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/R code/"
outPath       <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/"
plotPath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/Plots input/"
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/","/media/n/",path)
  inPath      <- sub("W:/","/media/n/",inPath)
  codePath    <- sub("W:/","/media/n/",codePath)
  outPath     <- sub("W:/","/media/n/",outPath)
}


source(paste(codePath,'functions.r',sep=""))

# load the true and observed stocks at the start of the simulation  from :
RecRegime <-  "Bayesian"

#define year ranges
ShortT    <-  ac(2013:2017)
MidT      <-  ac(2018:2027)
LongT     <-  ac(2027:2031)



##-------------------------------------------------------------------------------
## Setup array to save results
##-------------------------------------------------------------------------------



diags<-data.frame(scenario=NA,RecRegime=NA,Btrigger=NA,Blim=NA,Ftarget=NA,
      Risk2ShortT=NA,Risk2MidT=NA,Risk2LongT=NA,
      SSBend=NA,meanSSBLongTerm=NA,
      Fend=NA,meanF=NA,meanFLongTerm=NA,
      meanYieldShortTerm=NA,meanYieldMidTerm=NA,meanYieldLongTerm=NA,
      meanrelTACIAV=NA,noIAVrestrictup=NA,noIAVrestrictdown=NA,TACup=NA,TACdown=NA)                                                                                                  


##-------------------------------------------------------------------------------
## Load results
##-------------------------------------------------------------------------------


 counter <- 1
for (opt in 1:5)
{
scen          <- c("LTMP")         
TACvarlim     <- T                 
Fvarlim       <- T                 
BBscen        <- "AlternateBank"   
LastRecOp     <- "RCT3"            

sc.name<-paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep="")

outPath2<-paste(outPath,sc.name,"/",sep="")
source(paste(codePath,"07_scenarioDescription.r", sep=""))
mpPoints      <- get(scen)[[which(names(get(scen))==paste("opt",opt,sep=""))]]







load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalf.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalbiol.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalstocks.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalpercievedstocks.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalTAC.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalfSTF.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_Finalfishery.RData",           sep=""))
load(file=paste(outPath2,"/",scen,opt,mpPoints$FadultA,"_FinalmpPoints.RData",           sep=""))
source(paste(codePath,"functions.r",            sep=""))
source(paste(codePath,"04_forecastScenarios.r", sep=""))
load(file=paste(outPath,"settings.RData",       sep=""))
for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])



projPeriod           <- 2013:(futureMaxYr-1)  
print(counter)
print("question :")
print("is this better to use median or mean among stock replicates")
print("given that things may be multimodal due to different SR model")   
Ref<-mpPoints
diags[counter,"scenario"]   <- opt
diags[counter,"RecRegime"]   <- RecType
diags[counter,"Btrigger"]   <-Ref$Btrigger
diags[counter,"Blim"]   <-Ref$Blim
diags[counter,"Ftarget"]   <-Ref$Ftarget


##-------------------------------------------------------------------------------
## Diagnostics on results
##-------------------------------------------------------------------------------

Btrg  <- Ref$Btrigger
Blim  <- Ref$Blim
Bpa  <- Ref$Bpa


Ssb<-ssbb(biol,f,stockstore)
percSsb<-ssb(stockstore)
Fbar<-quantMeans((f[ac(4:8),]))
Fbar2<-quantMeans((harvest(stockstore)[ac(4:8),]))

# risk related to Blim
diags[counter,"Risk2ShortT"]             <- length(unique(which(Ssb[,ShortT]<Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100       # percentage of iteration that reach Blim
diags[counter,"Risk2MidT"]               <- length(unique(which(Ssb[,MidT]  <Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100
diags[counter,"Risk2LongT"]              <- length(unique(which(Ssb[,LongT] <Blim,arr.ind=T)[,"dim6"]))/dims(Ssb)$iter  *100


# stock and fishing mortality
diags[counter,"SSBend"]               <- round(median(c(apply(Ssb[,ac(futureMaxYr-1)],3:6,mean,na.rm=T))))   # or round(iterMeans(Ssb[,ac(futureMaxYr-1)]))
diags[counter,"meanSSBLongTerm"]      <- round(median(c(apply(Ssb[,LongT],3:6,mean,na.rm=T))))

diags[counter,"meanF"]                <- round(  median(  yearMeans(Fbar[,ac(projPeriod)])@.Data) ,3)
diags[counter,"meanFLongTerm"]        <- round( median(   yearMeans(Fbar[,LongT])@.Data) ,3)
diags[counter,"Fend"]                 <- round( median   (Fbar[,ac(futureMaxYr-1)]@.Data) ,3)

# difference between percieved and true stocks
diags[counter,"SSBabsBias"]           <-  round(apply( yearMeans(100*abs(percSsb[,ac(projPeriod)]-Ssb[,ac(projPeriod)])/Ssb[,ac(projPeriod)]),1:5,median,na.rm=T),3)
diags[counter,"SSBBias"]              <-  round(apply( yearMeans(100*(percSsb[,ac(projPeriod)]-Ssb[,ac(projPeriod)])/Ssb[,ac(projPeriod)]),1:5,median,na.rm=T),3)
diags[counter,"FbarabsBias"]          <-  round(apply( yearMeans(100*abs((Fbar2[,ac(projPeriod)])-(Fbar[,ac(projPeriod)]))/(Fbar[,ac(projPeriod)])),1:5,median,na.rm=T),3)
diags[counter,"FbarBias"]             <-  round(apply( yearMeans(100*((Fbar2[,ac(projPeriod)])-(Fbar[,ac(projPeriod)]))/(Fbar[,ac(projPeriod)])),1:5,median,na.rm=T),3)



# catches and quotas
diags[counter,"meanYieldShortTerm"]   <- round(median(c(yearMeans((computeLandings(fishery)[,ShortT])))))
diags[counter,"meanYieldMidTerm"]     <- round(median(c(yearMeans((computeLandings(fishery)[,MidT])))))
diags[counter,"meanYieldLongTerm"]    <- round(median(c(yearMeans((computeLandings(fishery)[,LongT])))))




# quota variability
diags[counter,"meanrelTACIAV"]        <- round(median(c(apply(abs(TAC[,ac(projPeriod[2]:rev(projPeriod)[1])] - TAC[,ac(projPeriod[1]:rev(projPeriod)[2])]) /TAC[,ac(projPeriod[2]:rev(projPeriod)[1])] * 100,3:6,mean,na.rm=T))),3)
IAVUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] == 1.2* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
IAVDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] == 0.8* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
#  #- Average number of times the IAV rule is applied upwards or downwards
diags[counter,"noIAVrestrictup"]<-0
if((nrow(IAVUp)) > 0 ){
  a <- IAVUp
  diags[counter,"noIAVrestrictup"]    <- max(0,median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x),na.rm=T)
}
 diags[counter,"noIAVrestrictdown"]<-0
if((nrow(IAVDown)) > 0 ){
  a <- IAVDown
  diags[counter,"noIAVrestrictdown"]  <- max(0,median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x),na.rm=T)
}
#
#  #- Which TAC of the runs go up and which go down
resUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
resDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
#
#  #- Mean increase in TAC is TAC goes up, or mean decrease in TAC is TAC goes down
diags[counter,"TACup"]   <- round(mean((TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] - TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])]))
diags[counter,"TACdown"] <- round(mean((TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])] - TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])]))
#




counter <- counter + 1
}


write.csv(t(diags),file=paste(outPath,"/tables_diags",paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep=""),".csv",sep=""),row.names=T)
