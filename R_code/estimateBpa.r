rm(list=ls())

library(FLSAM)
library(MASS)
library(msm)
wine <- F

path          <- "N:/Projecten/ICES WG/WKHELP/"
inPath        <- "N:/Projecten/ICES WG/WKHELP/data/"
codePath      <- "N:/Projecten/ICES WG/WKHELP/R/"
outPath       <- "N:/Projecten/ICES WG/WKHELP/Results/"
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("N:/","/media/n/",path)
  inPath      <- sub("N:/","/media/n/",inPath)
  codePath    <- sub("N:/","/media/n/",codePath)
  outPath     <- sub("N:/","/media/n/",outPath)
}
load(file=paste(outPath,"NSH.RData",            sep=""))


#- Settings
histMinYr   <- 1947
histMaxYr   <- 2011
futureMaxYr <- histMaxYr + 11
histPeriod  <- ac(histMinYr:histMaxYr)
projPeriod  <- ac((histMaxYr+1):futureMaxYr)
recrPeriod  <- ac(2003:2011)
selPeriod   <- ac(1997:2011)
fecYears    <- ac(2003:2011)
nyrs        <- futureMaxYr-dims(NSH)$maxyear
nits        <- 1000

source(paste(codePath,"functions.r",sep=""))

load(file.path("~/WKHELP"/,"iter_retro.RData"))
load(file.path(outPath,"catchSurveys.RData"))

iter_ssbFLQ <- FLQuant(array(NA,dim=c(1,length(histPeriod),1,1,nyrs,nits)),dimnames=list(age="all",year=histPeriod,unit="unique",season="all",area=1:nyrs,iter=1:nits))
iter_errorSSB <- iter_ssbFLQ[,,,,1:10,]
for(iRun in 1:nits){
  print(iRun)
  iter_ssb                <- lapply(iter_retro[[iRun]],function(x){return(window(ssb(x + NSH),start=histMinYr,end=histMaxYr))})
  iter_ssbFLQ[,,,,,iRun]  <- FLQuant(array(unlist(iter_ssb),dim=c(1,length(histPeriod),1,1,nyrs,1)),dimnames=list(age="all",year=histPeriod,unit="unique",season="all",area=1:nyrs,iter=iRun))
  for(i in 2:11)
    iter_ssbFLQ[,ac(histMaxYr-i+2),,,i,iRun] <- NA
  iter_error  <- exp(sweep(log(iter_ssbFLQ[,,,,2:11,iRun]),c(1:4,6),log(ssb(NSH)[,histPeriod,,,,1]),"-"))
  for(i in 1:10)
    iter_errorSSB[,(1+i):65,,,i,iRun] <- iter_error[,1:(65-i),,,i,]
}
    

iter_errorSSB2 <- FLQuant(array(NA,dim=c(1,length(histPeriod),1,1,1,nits*nyrs)),dimnames=list(age="all",year=histPeriod,unit="unique",season="all",area=1,iter=1:(nits*nyrs)))
counter <- 1
for(iRun in 1:1000){
  for(iRetro in 1:10){
    iter_errorSSB2[,,,,,counter] <- iter_errorSSB[,,,,iRetro,iRun]
    counter <- counter + 1
    print(counter)
  }
}
dat <- as.data.frame(iter_errorSSB2)
boxplot(dat$data ~ as.factor(dat$year),xlab="Years (aligned)",ylab="Fraction different from 'truth'")
abline(h=1,col=2,lty=2,lwd=2)
qts <- aggregate(dat$data,by=list(dat$year),quantile,probs=c(0.05,0.95),na.rm=T)
lines(y=qts$x[,1],x=as.factor(1947:2011),lty=2,lwd=2,col="blue")
lines(y=qts$x[,2],x=as.factor(1947:2011),lty=2,lwd=2,col="blue")

