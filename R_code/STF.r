#
# RUN with R2.13.2    (32-bit)
# packages :
# FLCore 2.4
# FLAssess 2.4
# FLSAM 0.99-991
#


rm(list=ls())
library(FLCore)
library(FLAssess)



work.dir    <-  "M:/WGWIDE/NEA Benchmark 2014/update advice"
output.dir  <-  "M:/WGWIDE/NEA Benchmark 2014/update advice/STFresults" 
setwd(work.dir) 
assess.name <-  "NEAMack-for-update-advice-2014"



source("flsam.r")
run.dir<-paste(work.dir,"/",assess.name,"/run",sep="")
source("SAM2FLR.r")

save(Mac.sam,file="mac-assessment.RData")

library(FLSAM)
load("mac-assessment.RData")                 # load the FLSAM object with tha assessment

# defines the time frame
stY <- dims(Mac.sam)$minyear
TaY <- dims(Mac.sam)$maxyear-1  #Terminal assessment year
ImY <- TaY+1                #Intermediate Year
AdY <- TaY+2                #Advice year
CtY <- TaY+3                #Continuation year - not of major concern but used in calculations in places
tbl.yrs     <- as.character(c(ImY,AdY,CtY))   #Years to report in the output table



# creates an FLStock based on assessment input and output data
Mac<-FLStock(stock.n(Mac.sam))
Mac<-Mac+Mac.sam
Mac@catch.n@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(work.dir,assess.name,"data","cn.dat"),skip=5))
Mac@catch.wt@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(work.dir,assess.name,"data","cw.dat"),skip=5))
Mac@catch<-computeCatch(Mac)
Mac@landings.n<-Mac@catch.n
Mac@landings.wt<-Mac@catch.wt
Mac@landings<-Mac@catch
Mac@stock.wt@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(work.dir,assess.name,"data","sw.dat"),skip=5))[,-dim(Mac@stock.wt)[2]]
m(Mac)<-0.15
Mac@mat@.Data[,ac(stY:TaY),,,,]<-t(read.table(file.path(work.dir,assess.name,"data","mo.dat"),skip=5))[,-dim(Mac@stock.wt)[2]]
Mac@harvest.spwn@.Data[,ac(stY:ImY),,,,]<-t(read.table(file.path(work.dir,assess.name,"data","pf.dat"),skip=5))#[,-dim(Mac@stock.wt)[2]]
Mac@m.spwn@.Data[,ac(stY:ImY),,,,]<-t(read.table(file.path(work.dir,assess.name,"data","pm.dat"),skip=5))
Mac@range[6]<-4
Mac@range[7]<-8



# prepare init file for RCT3
R<-stock.n(Mac)[ac(0),ac(1990:2012)]
R<-c(R@.Data)
R[length(R)]  <- -11
IBTS.index<-c(rep(-11,8),c(
0.467656705,
0.638432004,
0.232791963,
0.638170742,
0.65235593 ,
0.358648186,
0.809318184,
1.190326879,
1.063104797,
0.385761958,
0.636563318,
0.330024847,
0.503952216,
0.92553681 ,
0.834884058
))                         # this is the log transformed IBTS index
years<-1990:TaY


# remove files in the RCT3 folder !!!!

write.table(data.frame("RCT3 for NEA Mackerel"),file="RCT3/RCT3init.txt",quote=F,col.names=FALSE,row.names=FALSE,append=TRUE,sep="\t")
write.table(data.frame(1,length(R),2,"SAM","IBTS.index"),file="RCT3/RCT3init.txt",quote=F,col.names=FALSE,row.names=FALSE,append=TRUE,sep="\t")
write.table(data.frame(years,R,IBTS.index),file="RCT3/RCT3init.txt",col.names=FALSE,quote=F,row.names=FALSE,append=TRUE,sep="\t")
write.table(data.frame(c("SAM","IBTS.index")),file="RCT3/RCT3init.txt",col.names=FALSE,quote=F,row.names=FALSE,append=TRUE,sep="\t")

source("RCT3/RCT3v4a.r")
Rct3<-RCT3("RCT3/RCT3init.txt",logged=T)
Rct3$plot()
Rct3$input()
RCT3res<-Rct3$output()
Rct3$summary()
Rct3$printout("RCT3/RCT3output.txt")




# Deal with recruitment :
# terminal year recruitment to be replaced by output of RCT3
stock.n(Mac)[1,ac(TaY)]     <-      RCT3res$Years$WAPred  # RCT3 output using a predictor the IBTS index for the terminal assessmnet year
# recompute the survivors (year ImY) at age 1
stock.n(Mac)[2,ac(ImY)]    <- stock.n(Mac)[1,ac(TaY)] * exp (-(harvest(Mac)[1,ac(TaY)]+m(Mac)[1,ac(TaY)]))

# define the recruitment for the short term forecast (starting at ImY )
rec.years <- (1990:(TaY-1))
Mac.srr <- sr(as.FLSR(trim(Mac,year=rec.years),model=geomean))         # default is geomean used
# if desired, replace the geomean by the tapered time weighted mean from RCT 3
replace.by.tapered.time.weighted.mean <- F
tapered.time.weighted.mean       <- exp(RCT3res$Intercept[2])
if(replace.by.tapered.time.weighted.mean) params(Mac.srr)[1]  <- tapered.time.weighted.mean


##compute the survivors
#survivors<-stock.n(Mac)[,ac(ImY)]
#survivors[ac(1:12),]  <-   stock.n(Mac)[ac(0:11),ac(TaY)] * exp (-(harvest(Mac)[ac(0:11),ac(TaY)]+m(Mac)[ac(0:11),ac(TaY)]))
#survivors[ac(12),]    <-   survivors[ac(12),]  + stock.n(Mac)[ac(12),ac(TaY)] * exp (-(harvest(Mac)[ac(12),ac(TaY)]+m(Mac)[ac(12),ac(TaY)]))
#survivors[ac(0),]     <-   params(Mac.srr)[1] 
#
#
#ns<-stock.n(Mac)
#ns[]<-NA
#ns[ac(1:12),ac(1981:2013)]  <-   stock.n(Mac)[ac(0:11),ac(1980:2012)] * exp (-(harvest(Mac)[ac(0:11),ac(1980:2012)]+m(Mac)[ac(0:11),ac(1980:2012)]))
#ns[ac(12),ac(1981:2013)]    <-   ns[ac(12),ac(1980:2012)]  + stock.n(Mac)[ac(12),ac(1980:2012)] * exp (-(harvest(Mac)[ac(12),ac(1980:2012)]+m(Mac)[ac(12),ac(1980:2012)]))
#
#ns<-propagate(ns,2)
#ns[,,,,,2]<-stock.n(Mac)


#Expand stock object
#NEA.Mac.proj <- stf(NEA.Mac,nyears=4,wts.nyears=3,arith.mean=TRUE,na.rm=TRUE)
Mac.proj <- stf(trim(Mac, year=(TaY-5):TaY),nyears=3,wts.nyears=3,arith.mean=TRUE,na.rm=TRUE,fbar.nyears=3) # to use object saved before rounding

# use the actual 2013 values for F prop and M prop
harvest.spwn(Mac.proj)[,ac(ImY)]      <-  harvest.spwn(Mac)[,ac(ImY)]
harvest.spwn(Mac.proj)[,ac(AdY:CtY)]  <-  yearMeans(harvest.spwn(Mac.proj)[,ac((TaY-1):ImY)])
m.spwn(Mac.proj)[,ac(ImY)]    <-  m.spwn(Mac)[,ac(ImY)]
m.spwn(Mac.proj)[,ac(AdY:CtY)]  <-  yearMeans(m.spwn(Mac.proj)[,ac((TaY-1):ImY)])

# update the survivors in the stf object
Mac.proj@stock.n[,ac(ImY)]  <- stock.n(Mac)[,ac(ImY)]
Mac.proj@stock.n[1,ac(c(ImY,AdY,CtY))] <- params(Mac.srr)[1]



# do the forecast over the intermediate year first

#create a the function necessary
returnCatch<-function(Fmult,stf,yr)
{

}

Catch2F<-function(stf,Catch)
{



#Setup options to suite the advice sheet  the choice here is based on the following:
#0 catch,  role over + and - 20% on TAC, F=0.20,0.21,0.22 from management plan  
ImY.catch <- 895336  # estimated catches expected due to TAC, discards, payback, overfishing, unilateral quotas etc.

ImY.TAC <- 895336  


options.l <- list(#Zero catch
                  "Catch(2014) = Zero"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity="catch",
                                          val=c(ImY.catch,0,0))),
                  #2010 Catch is XXX, followed by -20% Catch reduction => XXX
                  "Catch(2014) = 2013 catch (excl. interannual transfer and discard) -20%"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","catch","f"),
                                          rel=c(NA,NA,AdY),
                                          val=c(ImY.catch,ImY.TAC*0.80,1))),
                  #2009 and 2010 Catch is XXX
                  "Catch(2014) = 2013 catch (excl. interannual transfer and discard)"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","catch","f"),
                                          rel=c(NA,NA,AdY),
                                          val=c(ImY.catch,ImY.TAC,1))),
                  #2010 Catch is XXX, followed by +20% Catch increase => 2011 Catch XXX
                  "Catch(2014) = 2013 catch (excl. interannual transfer and discard) +20%"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","catch","f"),
                                          rel=c(NA,NA,AdY),
                                          val=c(ImY.catch,ImY.TAC*1.20,1))),
                 #2010 Catch is XXX, followed Fbar= 0.20
                  "Fbar(2014) = 0.20"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","f","f"),
                                          val=c(ImY.catch,0.20,0.20))),
                 #2010 Catch is XXX, followed Fbar= 0.21
                  "Fbar(2014) = 0.21"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","f","f"),
                                          val=c(ImY.catch,0.21,0.21))),
                 #2010 Catch is XXX, followed Fbar= 0.22
                  "Fbar(2014) = 0.22 (Fmsy)"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","f","f"),
                                          val=c(ImY.catch,0.22,0.22))),
                 #
                  "Fbar(2014) = 0.23 (Fpa)"=
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","f","f"),
                                          val=c(ImY.catch,0.23 ,0.23 ))),
                                          
 "Fbar(2014) = F  msy "=                   # 0.4*F2010+0.6*FMSY ===  0.4*0.2772+0.6*0.23 ===  0.24888
                    projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","f","f"),
                                          val=c(ImY.catch,0.25 ,0.25 )))                                          
) #End options list

#Multi-options table - standard one to show wider range of options for the report
# F multipliers from 0 to 2 *roll over F
fmult.targs  <- seq(0,2,by=0.1)
mult.opts.l <- lapply(as.list(fmult.targs),function(fmult) {
                          projectControl(data.frame(year=c(ImY,AdY,CtY),
                                          quantity=c("catch","f","f"),
                                          rel=c(NA,ImY,AdY),
                                          val=c(ImY.catch,fmult,1)))
                  })
names(mult.opts.l) <- sprintf("Fmult(2012) = %4.3f",fmult.targs)

#Calculate options for two option tables
Mac.options   <- lapply(options.l,function(ctrl) {project(Mac.proj,ctrl,Mac.srr)})
names(Mac.options)<-names(options.l)
Mac.mult.opts <- lapply(mult.opts.l,function(ctrl) {project(Mac.proj,ctrl,Mac.srr)})
names(Mac.mult.opts)<-names(mult.opts.l)


### ======================================================================================================
### Write Options Tables
### ======================================================================================================

output.base <- output.dir
#Document input settings
input.tbl.file <-paste(output.base,"/options - input.csv",sep=".")
write.table(NULL,file=input.tbl.file,col.names=FALSE,row.names=FALSE)
input.tbl.list <- list(N="stock.n",M="m",Mat="mat",PF="harvest.spwn",
                       PM="m.spwn",SWt="stock.wt",Sel="harvest",CWt="catch.wt")
for(yr in c(ImY,AdY,CtY)){
    col.dat <- sapply(input.tbl.list,function(slt) slot(Mac.proj,slt)[,as.character(yr),drop=TRUE])
    write.table(yr,file=input.tbl.file,col.names=FALSE,row.names=FALSE,append=TRUE,sep=",")
    write.table(t(c("Age",colnames(col.dat))),file=input.tbl.file,col.names=FALSE,row.names=FALSE,append=TRUE,sep=",")
    write.table(col.dat,file=input.tbl.file,col.names=FALSE,row.names=TRUE,append=TRUE,sep=",",na="-")
    write.table("",file=input.tbl.file,col.names=FALSE,row.names=FALSE,append=TRUE,sep=",")
}

#Detailed options table
options.file <-paste(output.base,"/options - details.csv",sep=".")
write.table(NULL,file=options.file,col.names=FALSE,row.names=FALSE)
for(i in 1:length(Mac.options)) {
    opt <- names(Mac.options)[i]
    stk <- Mac.options[[opt]]
    #Now the F and N by age
    nums.by.age <- stk@stock.n[,tbl.yrs,drop=TRUE]
    colnames(nums.by.age) <- sprintf("N(%s)",tbl.yrs)
    f.by.age    <- stk@harvest[,tbl.yrs,drop=TRUE]
    colnames(f.by.age) <- sprintf("F(%s)",tbl.yrs)
    age.tbl     <- cbind(Age=rownames(f.by.age),N=nums.by.age,F=f.by.age)
    #And now the summary tbl
    sum.tbl     <- cbind(Year=tbl.yrs,SSB=ssb(stk)[,tbl.yrs],
                        F.bar=fbar(stk)[,tbl.yrs],Yield=computeCatch(stk)[,tbl.yrs])
    #Now, bind it all together
    sum.tbl.padding <- matrix("",nrow=nrow(age.tbl)-nrow(sum.tbl),ncol=ncol(sum.tbl))
    comb.tbl    <- cbind(age.tbl," ",rbind(sum.tbl,sum.tbl.padding))
    #And write it - hdr first, then the rest
    write.table(sprintf("%s). %s",letters[i],opt),options.file,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
    write.table(t(colnames(comb.tbl)),options.file,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
    write.table(comb.tbl,options.file,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
    write.table(c(""),options.file,append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
}

#Options summary table
opt.sum.tbl <- function(stcks,fname) {
options.sum.tbl <- sapply(as.list(1:length(stcks)),function(i) {
opt <- names(stcks)[i]
stk <- stcks[[opt]]
#Build up the summary
sum.tbl <- data.frame(Rationale=ac(names(stcks)[i]),
F.ImY=fbar(stk)[,as.character(ImY),drop=TRUE],
Catch.ImY=computeCatch(stk)[,as.character(ImY),drop=TRUE],
SSB.ImY=ssb(stk)[,as.character(ImY),drop=TRUE],

TSB.ImY=quantSums(stock.wt(stk)*stock.n(stk))[,as.character(ImY),drop=TRUE],

F.AdY=fbar(stk)[,as.character(AdY),drop=TRUE],
Catch.AdY=computeCatch(stk)[,as.character(AdY),drop=TRUE],
SSB.AdY=ssb(stk)[,as.character(AdY),drop=TRUE],

TSB.AdY=quantSums(stock.wt(stk)*stock.n(stk))[,as.character(AdY),drop=TRUE] ,

SSB.CtY=ssb(stk)[,as.character(CtY),drop=TRUE],
TSB.CtY=quantSums(stock.wt(stk)*stock.n(stk))[,as.character(CtY),drop=TRUE] )

})
options.sum.tbl <- t(options.sum.tbl)
options.sum.tbl[,1]<-names(stcks)
colnames(options.sum.tbl) <- c("Rationale",
sprintf("Fbar (%i)",ImY),sprintf("Catch (%i)",ImY),sprintf("SSB (%i)",ImY),sprintf("TSB (%i)",ImY),
sprintf("Fbar (%i)",AdY),sprintf("Catch (%i)",AdY),sprintf("SSB (%i)",AdY),sprintf("TSB (%i)",AdY),
sprintf("SSB (%i)",CtY),sprintf("TSB (%i)",CtY) )
write.csv(options.sum.tbl,file=fname,row.names=FALSE)
}
opt.sum.tbl(stcks=Mac.options,fname=paste(output.base,"/options - summary.csv",sep="."))
opt.sum.tbl(stcks=Mac.mult.opts,fname=paste(output.base,"/multi-options - summary.csv",sep="."))
















