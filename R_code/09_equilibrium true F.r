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

for (fequ in seq(0.0,0.6,0.01))
{              
              #- Load objects
              load(file=paste(outPath,"Mac.RData",            sep=""))
              load(file=paste(outPath,"Macctrl.RData",        sep=""))
              load(file=paste(outPath,"biol.RData",           sep=""))
              load(file=paste(outPath,"fishery.RData",        sep=""))
              load(file=paste(outPath,"propN.RData",        sep=""))
              load(file=paste(outPath,"propWt.RData",       sep=""))
              load(file=paste(outPath,"ctch.RData",         sep=""))
              load(file=paste(outPath,"surveys.RData",        sep=""))
              load(file=paste(outPath,"stocks.RData",         sep=""))
              load(file=paste(outPath,"settings.RData",       sep=""))
              load(file=paste(outPath,"resNFinal.RData",      sep=""))
              load(file=paste(outPath,"resFFinal.RData",      sep=""))
              load(file=paste(outPath,"SRmod.RData",      sep=""))
              load(file=paste(outPath,"resFFinal.RData",      sep=""))
              for(i in 1:length(settings)) assign(x=names(settings)[i],value=settings[[i]])
              
              
              source(paste(codePath,"functions.r",            sep=""))
              source(paste(codePath,"04_forecastScenarios.r", sep=""))
              
              RecType<-settings$RecType
              
              
              #   nits<-3
              
                #------------------------------------------------------------------------------#
                # 0) setup TACS & F's and Survivors and maximum change in effort
                #------------------------------------------------------------------------------#
              maxEff                                <- 1000
              TAC                                   <- FLQuant(NA,dimnames=list(age="all",year=histMinYr:(futureMaxYr+3),unit=c("A"),season="all",area=1,iter=1:nits))
              TAC[,ac(2000:2013),"A"]               <- c(612,670,683,583,532,422,444,502,458,605,885,959,927,906) *1000
              #TACusage                              <- FLQuant(array(rep(c(rep(1,length(histMinYr:futureMaxYr)+3),rep(0.539755,length(histMinYr:futureMaxYr)+3),
              #                                                   rep(1,length(histMinYr:futureMaxYr)+3),rep(1,length(histMinYr:futureMaxYr)+3)),nits),dim=c(1,length(histMinYr:futureMaxYr)+3,4,1,1,nits)),dimnames=dimnames(TAC[,ac(histMinYr:(futureMaxYr+3))]))
              TACusage                              <- FLQuant(array(rep(c(rep(1,length(histMinYr:futureMaxYr)+3)),nits),dim=c(1,length(histMinYr:futureMaxYr)+3,1,1,1,nits)),dimnames=dimnames(TAC[,ac(histMinYr:(futureMaxYr+3))]))
              
              HCRTAC                                <- TAC; HCRTAC[] <- NA; SSB <- HCRTAC[,,1]; HCRSSB <- SSB
              
              stockstore <- stocks # object to store the perceived stocks. The object "stocks" being updated at each time step, it does keep track of the percieved stock
              
              f                                     <- FLQuant(NA,dimnames=dimnames(fishery@landings.n))
              for(iFsh in dimnames(f)$unit)
                f[,ac(histMinYr:histMaxYr),iFsh]    <- sweep(harvest(stocks)[,ac(histMinYr:histMaxYr)],c(1,3:5),propN[,,iFsh],"*")
              fSTF                                  <- f; fSTF@.Data[] <- NA
              survivors                             <- FLQuant(NA,dimnames=dimnames(n(biol)[,ac(2011)]))
              
                #------------------------------------------------------------------------------#
                # 1) Select Scenario's   not relevant  here
                #------------------------------------------------------------------------------#
              
              outPath2       <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/equilibrium/true/"
              
              scen          <- "equilTrue"              # 
              opt           <- ""                      # for multiple scenario combinations a counter
              #TACvarlim     <- T                      # whether or not to apply the 20% limit on TAC change
              #Fvarlim       <- T                      # whether or not to apply the 10% limit on Fbar change
              #BBscen        <- "AlternateBank"        # banking borrowing options :
              #                                                     # "Banking"          : always bank 
              #                                                     # "Borrowing"        : always borrow
              #                                                     # "AlternateBank"    : bank first and then alternate
              #                                                     # "AlternateBorrow"  : borrow first and then alternate
              #                                                     # "MinVar"           : use BB to minise TAC variability
              #                                                     
              LastRecOp     <- "RCT3"                  # option to replace the last estimated recruitment : "SAM", "geom", "RCT3"
              #                                                     # "SAM"  = don't overwrite the SAM estimmate
              #                                                     # "geom" = replace by geomean 1990/(TaY-1)
              #                                                     # "RCT3" = replace by RCT3 output
              mpOptions<-list()
              #
              #sc.name<-paste(scen,opt,"_TACvarlim",TACvarlim,"_Fvarlim",Fvarlim,"_",BBscen,"_LastRec",LastRecOp,sep="")
              #
              #outPath2<-paste(outPath,sc.name,"/",sep="")
              source(paste(codePath,"07_scenarioDescription.r", sep=""))
              mpPoints      <- list(Fequ=fequ)
              #
              #dir.create(outPath2)
              
                #------------------------------------------------------------------------------#
                # 2) Start running
                #------------------------------------------------------------------------------#
              
              start.time <- Sys.time()
              for (iYr in an(projPeriod)){
                cat(iYr,"\n")
                cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
              
                #- Define mortality rates for iYr-1 to calculate survivors to iYr
                m           <- m(biol)[,ac(iYr-1),,]                                      
                z           <- unitSums(f[,ac(iYr-1),,,,]) + m
                ImY<-iYr
                #- Update biological model to iYr
                  #- Survivors
                survivors   <- n(biol)[,ac(iYr-1)] * exp(-z)
                n(biol)[ac((range(biol,"min")+1):range(biol,"max")),ac(iYr),,] <- survivors[-dim(survivors)[1],,,,,]@.Data
              
                  # - Recruitment
              ssbtp<-ssbb(biol[,ac(iYr-1),,,,],f[,ac(iYr-1),,,,],stockstore[,ac(iYr-1),,,,])
              for (its in 1:nits)   n(biol)[1,ac(iYr),,,,its][] <- B2R (ssbtp[,,,,,its],SRmod[its,],its,iYr,RecType)
              
            
              
                  #- Plusgroup
                if (!is.na(range(biol,"plusgroup"))){
                  n(biol)[ac(range(biol,"max")),ac(iYr),] <- n(biol)[ac(range(biol,"max")),ac(iYr),] + survivors[ac(range(biol,"max"))]
                }
              
                cat("\n Finished biology \n")
                cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
                   
                #- Update fishery to year iYr-1
                landings.n(fishery)[,ac(iYr-1)]     <- sweep(sweep(f[,ac(iYr-1),,,,],c(1:2,4:6),z,"/"),c(1:2,4:6),n(biol)[,ac(iYr-1)]*(1-exp(-z)),"*")
              
               
                TAC[,ac(iYr-1)] <-  quantSums(landings.wt(fishery)[,ac(iYr-1)] *  landings.n(fishery)[,ac(iYr-1)])
                
                
                cat("\n Finished forecast \n")
                cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
                
                #-Calculate effort accordingly (assuming constant catchability)
                f[,ac(ImY)]               <- sweep(catch.sel(fishery)[,ac(ImY)],c(3,6),mpPoints$Fequ,"*")
                cat("\n Finished effort calc \n")
                cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
              
                #- Save each year the output
                #save.image(file=paste(outPath,scen,"_",opt,"_",mpPoints$FadultA,"_",iYr,".RData",sep=""))
                #save.image(file=paste("/home/hintz001/WKHELP_test2_",iYr,".RData",sep=""))
                #save(file=paste("D:/WKHELP_test3_",iYr,".RData",sep=""))
              }
              
              stockstore@landings.n   <- stockstore@harvest * stockstore@stock.n * (1- exp(-stockstore@harvest - stockstore@m)) / (stockstore@harvest + stockstore@m)
              stockstore@landings.wt<-stockstore@catch.wt
                #-------------------------------------------------------------------------------
                # 3): Save the objects
                #-------------------------------------------------------------------------------
              



save(biol          ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalbiol.RData",        sep=""))
save(fishery       ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalfishery.RData",     sep=""))
save(f             ,file=paste(outPath2,scen,opt,ac(mpPoints$Fequ),"_Finalf.RData",           sep=""))

}


