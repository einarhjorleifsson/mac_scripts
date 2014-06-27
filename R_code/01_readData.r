#-------------------------------------------------------------------------------
# Mackerel management plan evaluation
#
# Author: Thomas Brunel    (based on Niels Hintzen code for North Sea herring management plan evaluation)
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
require(FLSAM)

if(substr(R.Version()$os,1,3) == "lin") 
  {
   path_prefix <- "~/ass/2014/mac/MSE"
   } else {
     path_prefix <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE"
   }

path          <- path_prefix
inPath        <- paste(path_prefix,"Data",sep="/")
codePath      <- paste(path_prefix,"R_code",sep="/")
outPath       <- paste(path_prefix,"Results",sep="/")

#if(substr(R.Version()$os,1,3)== "lin"){
#  path        <- sub("W:/","/media/n/",path)
#  inPath      <- sub("W:/","/media/n/",inPath)
#  codePath    <- sub("W:/","/media/n/",codePath)
#  outPath     <- sub("W:/","/media/n/",outPath)
#}

  #-------------------------------------------------------------------------------
  # 0): Read the data & assessment results (from the WKPELA 2014 benchmark assessment updated in April 2014 for the update advice)
  #-------------------------------------------------------------------------------

#- Read SAM output
setwd(path)

assess.name <-  "NEAMack-for-update-advice-2014"
# EINAR source(paste(codePath,"/utils/flsam.r",sep=""))
run.dir <- paste(inPath,"/",assess.name,"/run",sep="")
source(paste(codePath,"/utils/SAM2FLR.r",sep=""))
#library(FLSAM)
name(Mac.sam) <- "NEA Mackerel"


# defines the time frame
stY <- dims(Mac.sam)$minyear
TaY <- dims(Mac.sam)$maxyear-1  #Terminal assessment year
ImY <- TaY+1                #Intermediate Year

# creates an FLStock based on assessment input and output data
Mac<-FLStock(stock.n(Mac.sam))
Mac<-Mac+Mac.sam
Mac@catch.n@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(inPath,"/",assess.name,"data","cn.dat"),skip=5))
Mac@catch.wt@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(inPath,"/",assess.name,"data","cw.dat"),skip=5))
Mac@catch<-computeCatch(Mac)
Mac@landings.n<-Mac@catch.n
Mac@landings.wt<-Mac@catch.wt
Mac@landings<-Mac@catch
Mac@stock.wt@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(inPath,"/",assess.name,"data","sw.dat"),skip=5))[,-dim(Mac@stock.wt)[2]]
m(Mac)<-0.15
Mac@mat@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(inPath,"/",assess.name,"data","mo.dat"),skip=5))[,-dim(Mac@stock.wt)[2]]
Mac@harvest.spwn@.Data[,ac(stY:ImY),,,,]<-t(read.table(file.path(inPath,"/",assess.name,"data","pf.dat"),skip=5))#[,-dim(Mac@stock.wt)[2]]
Mac@m.spwn@.Data[,ac(stY:ImY),,,,]<-t(read.table(file.path(inPath,"/",assess.name,"data","pm.dat"),skip=5))
Mac@range[6]<-4
Mac@range[7]<-8
Mac@stock     <- computeStock(Mac)


#  survey indices (but not tags), this is just to have the recruitment index ready for the simulations
Mac.tun   <- readFLIndices(file.path(inPath,"/",assess.name,"data","survey.dat"))
Mac.tun[[1]]@catch.n[Mac.tun[[1]]@catch.n==-1]<-NA
Mac.tun[[3]]@catch.n[Mac.tun[[3]]@catch.n==-1]<-NA
Mac.tun[[1]]@type<-"biomass"
Mac.tun[[2]]@type<-"number"
Mac.tun[[3]]@type<-"number"


# control object
source(file.path(codePath,"setupControlObject.r"))

  #-------------------------------------------------------------------------------
  # 1): Save 2013 data
  #-------------------------------------------------------------------------------

save(Mac,       file=paste(outPath,"/Mac.RData",     sep=""))
save(Mac.sam,   file=paste(outPath,"/Macsam.RData",  sep=""))
save(Mac.tun,   file=paste(outPath,"/Mactun.RData",  sep=""))
save(Mac.ctrl,  file=paste(outPath,"/Macctrl.RData", sep=""))
