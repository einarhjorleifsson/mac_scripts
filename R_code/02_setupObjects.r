#-------------------------------------------------------------------------------
# Mackerel management plan evaluation
#
# Author: Thomas Brunel    (based on Niels Hintzen's code for North Sea herring management plan evaluation)
#         IMARES, The Netherland
#
# Performs an MSE of NEA Mackerel under different TAC scenario's
#
# Date: May/June-2014
#
# Build for R2.13.2, 32bits
# RUN with R2.13.2    (32-bit)
# packages :
# FLCore 2.4
# FLAssess 2.4
# FLSAM 0.99-991
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
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",path)
  inPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",inPath)
  codePath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
  outPath     <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",outPath)
    path        <- sub("WKPELA 2014 mackerel benchmark","WKPELA\ 2014\ mackerel\ benchmark",path)
  inPath      <- sub("WKPELA 2014 mackerel benchmark","WKPELA\ 2014\ mackerel\ benchmark",inPath)
  codePath    <- sub("WKPELA 2014 mackerel benchmark","WKPELA\ 2014\ mackerel\ benchmark",codePath)
  outPath     <- sub("WKPELA 2014 mackerel benchmark","WKPELA\ 2014\ mackerel\ benchmark",outPath)
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
load(file=paste(outPath,"Mac.RData",     sep=""))
load(file=paste(outPath,"Macsam.RData",  sep=""))
load(file=paste(outPath,"Mactun.RData",  sep=""))
load(file=paste(outPath,"/Macctrl.RData", sep=""))

#- Settings
assess.name <-  "NEAMack-for-update-advice-2014"
histMinYr   <- 1980
histMaxYr   <- 2012
nyrs        <- 40
futureMaxYr <- histMaxYr + nyrs
histPeriod  <- ac(histMinYr:histMaxYr)
projPeriod  <- ac((histMaxYr+1):futureMaxYr)
recrPeriod  <- ac(1990:2011)
selPeriod   <- ac(1990:2012)
fecYears    <- ac(1990:2011)
nits        <- 50
                    
                    
RecType     <-  "Bayesian"                # chose from
BiolType    <-  "ARMA"


settings    <- list(histMinYr=histMinYr,histMaxYr=histMaxYr,futureMaxYr=futureMaxYr,
                    histPeriod=histPeriod,projPeriod=projPeriod,recrPeriod=recrPeriod,
                    nyrs=nyrs,nits=nits,fecYears=fecYears,RecType=RecType,BiolType=BiolType)


source(paste(codePath,"functions.r",sep=""))


  #-------------------------------------------------------------------------------
  # 1): Create stock object & use vcov for new realisations
  #-------------------------------------------------------------------------------
Mac.sam@control<-Mac.ctrl
run.dir<-paste(inPath,"/",assess.name,"/run",sep="")

stocks                            <- monteCarloStock(Mac,Mac.sam,nits,run.dir=run.dir)
stocks                            <- window(stocks,start=histMinYr,end=futureMaxYr)
stocks@catch.n                    <- stocks@stock.n * stocks@harvest / (stocks@harvest + stocks@m) * (1 - exp(-stocks@harvest - stocks@m))
stocks@landings.n                 <- stocks@catch.n

stocks@harvest.spwn[,projPeriod]  <- stocks@harvest.spwn[,ac(histMaxYr)]
stocks@m.spwn[,projPeriod]        <- stocks@m.spwn[,ac(histMaxYr)]





  #-------------------------------------------------------------------------------
  # 2): create simulated data for the biology of the stock
  #-------------------------------------------------------------------------------

# a simple mean for sensitivity tests
if( BiolType == "MEAN")
{ 
stocks@mat      [,projPeriod][]               <-     yearMeans(stocks@mat[,ac((histMaxYr-2):histMaxYr)])                 ; print("no var in maturity")
stocks@stock.wt [,projPeriod][]               <-     yearMeans(stocks@stock.wt[,ac((histMaxYr-2):histMaxYr)])            ; print("no var in stock weights")
stocks@m        [,projPeriod][]               <-     0.15                                                               
}


# a simple mean for sensitivity tests
if( BiolType == "ARMA")
{ 

# load the objects previously generated
load(paste(inPath," ARMAfprop.RData",sep=""))
load(paste(inPath," ARMAmprop.RData",sep=""))
load(paste(inPath," ARMAstock.wt.RData",sep=""))
load(paste(inPath," ARMAcatch.wt.RData",sep=""))
load(paste(inPath," ARMAmat.RData",sep=""))


# adapt the size of the FLquants previously generated to the setup of this simulation
wt<-window(wt,start=an(projPeriod[1]),end=an(rev(projPeriod)[1]))
wt<-iter(wt,1:nits)
cw<-window(cw,start=an(projPeriod[1]),end=an(rev(projPeriod)[1]))
cw<-iter(cw,1:nits)
mat<-window(mat,start=an(projPeriod[1]),end=an(rev(projPeriod)[1]))
mat<-iter(mat,1:nits)
mprop<-window(mprop,start=an(projPeriod[1]),end=an(rev(projPeriod)[1]))
mprop<-iter(mprop,1:nits)
fprop<-window(fprop,start=an(projPeriod[1]),end=an(rev(projPeriod)[1]))
fprop<-iter(fprop,1:nits)

# overwritten all slots of the stocks objects
stocks@mat      [,projPeriod]               <-     mat
stocks@stock.wt [,projPeriod]               <-     wt
stocks@catch.wt [,projPeriod]               <-     cw
stocks@harvest.spwn [,projPeriod]           <-     fprop
stocks@m.spwn [,projPeriod]                 <-     mprop



stocks@m        [,projPeriod]               <-     0.15                                                               
}


 # quick check
                                  # apply(stocks@stock.wt,1:5,function(x) {sum(is.na(x))})
                                  # apply(stocks@catch.wt,1:5,function(x) {sum(is.na(x))})
                                  # apply(stocks@mat,1:5,function(x) {sum(is.na(x))})
                                  # apply(stocks@harvest.spwn,1:5,function(x) {sum(is.na(x))})
                                  # apply(stocks@m.spwn,1:5,function(x) {sum(is.na(x))})
                                  


  #-------------------------------------------------------------------------------
  # 3): Create survey object & use vcov for new realisations + error on realisations
  #-------------------------------------------------------------------------------
Rindex           <- lapply(Mac.tun,propagate,iter=nits)
Rindex           <-Rindex[[2]]

Rindex <- window(Rindex,start=range(Rindex)["minyear"],end=futureMaxYr)
dmns              <- dimnames(Rindex@index)
surv              <- FLQuant(NA,dimnames=dmns)

  #- Get redrawn survey Qs and Ks
load(file=file.path(run.dir, "random.param.RData"))

  #- Get the index of each parameter in the random.param object
Qidx              <- unlist(apply(Mac.ctrl@catchabilities,1,function(x)c(na.omit(x))))

  #- Create objects for surveyQ and surveyK's
surveyQ           <- FLQuants("R-idx(log transf)"=  FLQuant(NA,dimnames=dimnames(Rindex@index)))
surveyK           <- surveyQ
surveyK[[1]][]       <- 1


  #- Fill the Qs by survey
for(iYr in dimnames(surveyQ[[1]])$year)
surveyQ[["R-idx(log transf)"]][,iYr]      <- exp(random.param[,which(colnames(random.param) %in% "logFpar")[Qidx[grep("R-idx",names(Qidx))]]])

Rindex@index.q    <-  surveyQ[["R-idx(log transf)"]]  

  #- Index var no longer used but filled anyway
obsvar<-  exp(random.param[,which(colnames(random.param) %in% "logSdLogObs")[1+Qidx[grep("R-idx",names(Qidx))]]])
for(iYr in dimnames(surveyQ[[1]])$year)
Rindex@index.var[,iYr] <- obsvar
 
 
 
 
 
   #-------------------------------------------------------------------------------
  #- 4): catch estimation errors (assumed to be of the same magnitude as the residuas from the assessment for the period after 2000
  #-------------------------------------------------------------------------------

dmns              <- dimnames(trim(Mac@catch.n,year=1980:2012))
dmns$year         <- dmns$year[1]:futureMaxYr
dmns$iter         <- 1:nits
ctch              <- FLQuant(NA,dimnames=dmns)

  #- Take blocks of residuals, sample blocks from 1-10 and add up till length of timeseries
yrs       <- range(Mac)["minyear"]:futureMaxYr
saveBlcks <- matrix(NA,nrow=nits,ncol=length(yrs))
sam   <- sample(1:6,nits*length(yrs),replace=T)
samM  <- matrix(sam,nrow=nits,ncol=length(yrs))
cu<-apply(samM,1,cumsum)
cu[cu>length(yrs)] <-NA
dif<-(length(yrs)-cu)
bl<-samM*t(!is.na(cu))
bl[bl==0]<-NA

for (i in 1:nits)  
{
lst<-order(dif[,i])[1]
bl[i,lst]<-dif[lst,i]+bl[i,lst]

}
saveBlcks<-bl


  #- Take the sampled blocks and assign years to it
saveBlcksYrs <- array(NA,dim=c(nits,length(yrs),2),dimnames=list(nits=1:nits,yrs=yrs,strtstp=c("start","stop")))
for(iCol in 1:ncol(saveBlcksYrs)){
#  strstp  <- as.integer(runif(nits,range(Mac)["minyear"]+saveBlcks[,iCol]-1,range(Mac)["maxyear"]-saveBlcks[,iCol]+2))
  strstp  <- as.integer(runif(nits,2000+saveBlcks[,iCol]-1,(range(Mac)["maxyear"]-1)-saveBlcks[,iCol]+2))    # change the original code because we use only residuals since 2000
  rv      <- sample(c(T,F),nits,replace=T)
  strt    <- ifelse(rv==F,strstp,strstp-saveBlcks[,iCol]+1)
  stp     <- strt + saveBlcks[,iCol] - 1
  saveBlcksYrs[,iCol,"start"]  <- ifelse(rv==F,strt,stp)
  saveBlcksYrs[,iCol,"stop"]   <- ifelse(rv==F,stp,strt)
}

  #- Substract and calculate residuals (non-standardized)
Resids                                                        <- subset(residuals(Mac.sam),fleet=="Fleet 1")
iResids                                                       <- FLQuant(NA,dimnames=c(dimnames(stocks@stock.n)[1:5],iter="1"))

for(i in 1:nrow(Resids))
  iResids[ac(Resids$age[i]),ac(Resids$year[i]),]              <- exp(Resids$log.obs[i] - Resids$log.mdl[i])

  #- Fill the object with the residuals
for(iTer in 1:nits){
  blk <- which(is.na(saveBlcksYrs[iTer,,1])==F)
  idx <- ac(unlist(mapply(seq,from=saveBlcksYrs[iTer,blk,"start"],to=saveBlcksYrs[iTer,blk,"stop"])))
    #- Fill survey pattern with random draws of historic years
  iter(ctch[,      ac(yrs),],iTer)   <- iResids[,idx]
}
  #- Because residuals are not estimated everywhere, some are NA, replace with 1
ctch@.Data[which(is.na(ctch))] <- 1 
  
  
  
  #-------------------------------------------------------------------------------
  # 5): Create biological population object and define the SR models to be used
  #-------------------------------------------------------------------------------

biol                      <- as.FLBiol(stocks)


if(RecType=="geomean")
{
  #- Random draw from lognormal distribution for new recruitment, estimate lognormal parameters first
recrAge                   <- dimnames(rec(stocks))$age
pars                      <- optim(par=c(17.1,0.20),fn=optimRecDistri,recs=sort(c(rec(Mac[,ac(recrPeriod)]))),
                                  method="Nelder-Mead")$par
biol@n[1,projPeriod]      <- rtlnorm(length(projPeriod)*nits,mean=pars[1],sd=pars[2],lower=0.01*min(biol@n[recrAge,],na.rm=T))
}

if(RecType=="Bayesian") 
{
SRmod<-read.csv(file=paste(inPath,"SRbayes",recrPeriod[1],"_",rev(recrPeriod)[1],".csv",sep=""))[,-1]
cuales<- sample(1:1000,nits,replace=F)
SRmod<-SRmod[cuales,]
}

if(RecType=="EqSym")
{
SRmod<-read.csv(file=paste(inPath,"SREqSym",recrPeriod[1],"_",rev(recrPeriod)[1],".csv",sep=""))[,-1]
cuales<- sample(1:1000,nits,replace=F)
SRmod<-SRmod[cuales,]
}


if (RecType=="Indivit")
{
source(paste(codePath,"/LL for SR estimation.r",sep=""))
SRmod<-data.frame(mod=NA,A=NA,B=NA,sigma=NA)


updte=T
if(updte)  
{   
  for (iTer in 1:nits)
      {
      # stock recruit pairs for the corresponding iter
      R      <-  c(rec(stocks)[,ac(recrPeriod),,,,iTer]@.Data  )
      SSB    <-  c(ssb(stocks)[,ac(recrPeriod),,,,iTer]@.Data  )  
      source(paste(codePath,"SRfit function.r",sep=""))
      SRmod<-rbind(SRmod,data.frame(mod=res$model,A=res$a,B=res$b,sigma=res$sigma))
      cat(iTer,"\n") 
      write.csv(SRmod,file=paste(inPath,"SRIndivFit.csv",sep=""))
      }
  SRmod<-SRmod[-1,]    
  write.csv(SRmod,file=paste(inPath,"SRIndivFit.csv",sep=""))
  }



SRmod<-read.csv(file=paste(inPath,"SRIndivFit.csv",sep=""))
cuales<- sample(1:1000,nits,replace=F)
SRmod<-SRmod[cuales,]




 SRmod<-read.csv(file=paste(inPath,"SRIndivFit.csv",sep=""))
  SRmod<-SRmod[-1,]    
SRmod$mod<-factor(SRmod$mod)
SRmod<-SRmod[-1,]  
par(mfrow=c(3,2))
for (m in levels(SRmod$mod))
{
plot(density(SRmod$A[SRmod$mod==m]),main=paste(m,"A"))
plot(density(SRmod$B[SRmod$mod==m]),main=paste(m,"B"))
}



 aggregate(rep(1,dim(SRmod)[1]),list(SRmod$mod),sum)



}




  #-------------------------------------------------------------------------------
  # 6): compute the deviation from the "true" stock (Mac) for each replicate
  #-------------------------------------------------------------------------------
  
devN<- log(sweep(stock.n(stocks)[,ac(histPeriod)],1:5,stock.n(Mac)[,-34]  ,"/"))
devF<- log(sweep(harvest(stocks)[,ac(histPeriod)],1:5,harvest(Mac)[,-34]  ,"/"))

cvN<-(iterVars(devN))^0.5
cvF<-(iterVars(devF))^0.5



  
  # a): compute the autocorrelated part of the deviation 
  #-------------------------------------------------------------------------------
rhoN<-0.5  # assume an autocorrelation for the moment
rhoF<-0.5
#create an array which repeats the standard FLquant as many times as there are projection years

dnms<-dimnames(stock.n(stocks))
dnms$unit<-ac(c(2012,projPeriod))
En<-FLQuant(NA,dimnames=dnms)
Ef<-En

En[]<-           (1-rhoN^2)^0.5*rnorm(length(c(En[]@.Data)),0,1)
Ef[]<-           (1-rhoF^2)^0.5*rnorm(length(c(Ef[]@.Data)),0,1)

for (pYr in an(projPeriod))  for (a in 1:12) for (yr in  (histMinYr+1):futureMaxYr)  
                                                        {
                                                        En[ac(a),ac(yr),ac(pYr)] <-   rhoN* En[ac(a-1),ac(yr-1),ac(pYr-1)] + En[ac(a),ac(yr),ac(pYr)]
                                                        Ef[ac(a),ac(yr),ac(pYr)] <-   rhoF* Ef[ac(a-1),ac(yr-1),ac(pYr-1)] + Ef[ac(a),ac(yr),ac(pYr)]
                                                        cat(a,"/",yr,"\n",sep="")
                                                        }

#
#expl<-c()
#styr<-2020 
#for (a in 0:12) expl<-c(expl,En[ac(a),ac(styr+a),ac(styr+a),,,1])
#plot(expl)
#acf(expl)
#

  # b): add the autocorrelated part and the random deviation
  #-------------------------------------------------------------------------------

dnms<-dimnames(stock.n(stocks))
dnms$unit<-ac(projPeriod)
devN<-FLQuant(NA,dimnames=dnms)
devN[]<-1 #rnorm(length(devN[]),0,1)
devF<-devN

for (pYr in an(projPeriod))  
      for (a in 0:12) 
            for (yr in  (histMinYr:histMaxYr)+(pYr-2013))
                        {
                        devN[ac(a),ac(yr),ac(pYr),,,] <-   sweep(devN[ac(a),ac(yr),ac(pYr),,,],1:5,cvN[ac(a),ac(yr-(pYr-2013))],"*")
                        devF[ac(a),ac(yr),ac(pYr),,,] <-   sweep(devF[ac(a),ac(yr),ac(pYr),,,],1:5,cvF[ac(a),ac(yr-(pYr-2013))],"*")
                          cat(pYr,"/",a,"/",yr,"\n",sep="")
                        }


devN<-devN*En[,,-1,,,]
devF<-devF*Ef[,,-1,,,]
#
#rm(expl)
#expl<-c()
#styr<-2015
#for (a in 0:12) expl<-c(expl,devN[ac(a),ac(styr+a),ac(styr+a),,,1])
#plot(expl)
#acf(expl)
#                        






  #-------------------------------------------------------------------------------
  # 7): Create fisheries object
  #-------------------------------------------------------------------------------

dmns                      <- dimnames(m(biol))
dmns$unit                 <- c("A")
fishery                   <- FLCatch(price=FLQuant(NA,dimnames=dmns))
name(fishery)             <- "catches"
desc(fishery)             <- "NEA Mackerel"
fishery@range             <- range(biol)


    #-------------------------------------------------------------------------------
    #- Partial Ns per fleet and plusgroup setting : not use here for mackerel
    #-------------------------------------------------------------------------------
  dmns$year               <- ac(histMaxYr+1); dmns$iter <- 1;
  propN                   <- FLQuant(1,dimnames=dmns); propWt <- FLQuant(1,dimnames=dmns)
  propWt  <-propN
 
  #-Take single fleet weights and numbers and multiply by the proportions
for(iFsh in dimnames(fishery@landings.n)$unit){
  fishery@landings.n[,  ac(histMinYr:histMaxYr),iFsh]         <- Mac@landings.n[,  ac(histMinYr:histMaxYr)]
  fishery@landings.wt[, ac(histMinYr:histMaxYr),iFsh]         <- Mac@landings.wt[,  ac(histMinYr:histMaxYr)]
  fishery@discards.n[,  ac(histMinYr:histMaxYr),iFsh]         <- 0
  fishery@discards.wt[, ac(histMinYr:histMaxYr),iFsh]         <- 0
}
fishery@landings.n@.Data[is.infinite(fishery@landings.n)==T]  <- 0
fishery@landings.n@.Data[is.na(fishery@landings.n)==T]        <- 0
fishery@landings[,    ac(histMinYr:histMaxYr)]                <- computeLandings(fishery[,ac(histMinYr:histMaxYr)])
fishery@discards[,    ac(histMinYr:histMaxYr)]                <- computeDiscards(fishery[,ac(histMinYr:histMaxYr)])
                          #check: computeLandings(Mac) / window(unitSums(fishery@landings),1980,2012) #must equal 1
  # overwrite by the catch weights simulated
fishery@landings.wt                      <- stocks@catch.wt

  #-Calculate deterministic landing.sel
units(harvest(stocks))="f"
landings.sel(fishery)[,ac(histMinYr:histMaxYr)]        <- FLQuant(sweep(harvest(stocks[,ac(histMinYr:histMaxYr)]),2:6,
                                                                               fbar(stocks[,ac(histMinYr:histMaxYr)]),"/"),
                                                                         dimnames=dimnames(stocks[,ac(histMinYr:histMaxYr)]@stock.n))

catch.q(     fishery)[]                                       <- 1
discards.sel(fishery)[]                                       <- 0
fishery@discards.wt[]                                         <- 0
fishery@discards.n[]                                          <- 0

 
  #- Vary selectivity of fleet (add random walk to last year selection pattern)
  #-------------------------------------------------------------------------------
                
dmns                                                          <- dimnames(Mac@harvest)
dmns$year                                                     <- projPeriod
dmns$iter                                                     <- 1:nits
ages                                                          <- dimnames(stocks@stock.n)$age

#- Create random walk over Fs (as multiplier from last years selection pattern)
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nyrs,rep(0,length(ages)),cov(apply(log(Mac@harvest[,selPeriod,drop=T]),1,diff)))),
                                                                               dim=c(length(ages),nyrs,1,1,1,nits)),
                                                                         dimnames=dmns)
qtil                                                          <- quantile(c(wF),probs=c(0.025,0.975))
wF@.Data[which(wF<qtil[1] | wF>qtil[2])][]                    <- 0
rwF                                                           <- FLQuant(aperm(apply(wF,c(1,3:6),cumsum),c(2,1,3:6)),dimnames=dmns)
rwF                                                           <- sweep(rwF,c(1,3:5),apply(log(Mac@harvest[,selPeriod]),1,mean),"+")
fbarages                                                      <- ac(range(Mac)["minfbar"]:range(Mac)["maxfbar"])
landsel                                                       <- sweep(exp(rwF),c(2:6),apply(exp(rwF)[fbarages,],2:6,mean),"/")

for(iFsh in 1:dims(fishery)$unit)
  landings.sel(fishery)[,projPeriod,iFsh]                     <- FLQuant(sweep(landsel,c(1),propN[,,iFsh],"*"),
                                                                        dimnames=dimnames(stocks[,ac(projPeriod)]@stock.n))
 




  #-------------------------------------------------------------------------------
  # 4): Save the objects
  #-------------------------------------------------------------------------------
  
outPath<-paste(outPath,"50its40years",sep="")  
  
  save(biol          ,file=paste(outPath,"biol.RData",          sep=""))
  save(pars          ,file=paste(outPath,"recPars.RData",       sep=""))
  save(fishery       ,file=paste(outPath,"fishery.RData",       sep=""))
  save(propN         ,file=paste(outPath,"propN.RData",         sep=""))
  save(propWt        ,file=paste(outPath,"propWt.RData",        sep=""))
  save(ctch          ,file=paste(outPath,"ctch.RData",          sep=""))
  save(landsel       ,file=paste(outPath,"landsel.RData",       sep=""))
  save(Rindex        ,file=paste(outPath,"surveys.RData",       sep=""))
  save(stocks        ,file=paste(outPath,"stocks.RData",        sep=""))
  save(settings      ,file=paste(outPath,"settings.RData",      sep=""))
  save(CV_stock.n    ,file=file.path(outPath,"CVstockn.RData"))
  save(CV_harvest    ,file=file.path(outPath,"CVharvest.RData"))
  save(SRmod         ,file=file.path(outPath,"SRmod.RData"))
  save(devN          ,file=paste(outPath,"resNFinal.RData",     sep=""))
  save(devF          ,file=paste(outPath,"resFFinal.RData",     sep=""))
  save.image(         file=paste(outPath,"setup14092012.RData", sep=""))

