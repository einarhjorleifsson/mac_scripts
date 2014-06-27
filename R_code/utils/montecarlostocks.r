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


#- Load objects
load(file=paste(outPath,"Mac.RData",     sep=""))
load(file=paste(outPath,"Macsam.RData",  sep=""))
load(file=paste(outPath,"Mactun.RData",  sep=""))
load(file=paste(outPath,"/Macctrl.RData", sep=""))

#- Settings
assess.name <-  "NEAMack-for-update-advice-2014"
histMinYr   <- 1980
histMaxYr   <- 2012
futureMaxYr <- histMaxYr + 11
histPeriod  <- ac(histMinYr:histMaxYr)
projPeriod  <- ac((histMaxYr+1):futureMaxYr)
recrPeriod  <- ac(1980:2011)
selPeriod   <- ac(1990:2012)
fecYears    <- ac(1990:2011)
nyrs        <- futureMaxYr-dims(Mac)$maxyear
nits        <- 100
settings    <- list(histMinYr=histMinYr,histMaxYr=histMaxYr,futureMaxYr=futureMaxYr,
                    histPeriod=histPeriod,projPeriod=projPeriod,recrPeriod=recrPeriod,
                    nyrs=nyrs,nits=nits,fecYears=fecYears)

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


