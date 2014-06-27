# 


rm(list=ls())
Path<-"M:/WGWIDE/NEA Benchmark 2014/update advice"

setwd(Path)
library(FLCore)

source("flsam.r")
run.dir<-"M:/WGWIDE/NEA Benchmark 2014/update advice/NEAMack-for-update-advice-2014/run"
source("SAM2FLR.r")

save(Mac.sam,file="mac-assessment.RData")
