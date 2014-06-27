#-------------------------------------------------------------------------------
# Mackerel management plan evaluation
#
# Author: Thomas Brunel    (based on Niels Hintzen code for North Sea herring management plan evaluation)
#         IMARES, The Netherlands
#
# Performs an MSE of NEA Mackerel under different TAC scenario's  :
#   
# This script simulates biological parameteres (stocks and catch weigth, maturity, F and Mprop) using a ARMA model where predictions
# are forced to conect with the last years of the observations
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



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   define location, load objects, define the dimentions
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  !!!  requires R3.0.2 to run the  fArma library
rm(list=ls())
library(FLCore)
library(fArma)

inPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/exploration/"
assess.name <-  "NEAMack-for-update-advice-2014"
source(paste(codePath,'armaSearch.r',sep=""))


load("W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/stocks.RData")
load("W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Results/fishery.RData")

 
 histMinYr<-an(dimnames(stock.wt(stocks))$year[1])
 histMaxYr<-2012
 histPeriod<-  histMinYr:histMaxYr
 futureMaxYr<-2060
 projPeriod<-(histMaxYr+1):futureMaxYr
 nits<-1000
 ages<-ac(0:12)
 



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of stock weights 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 sw<-read.table(file.path(inPath,"/",assess.name,"data","sw.dat"),skip=5)
 sw$Year<-1980:2013
 sw<-sw[sw$Year!=2013,]
 
 w<-expand.grid(Year=1980:2012, Age=0:12)
 w$wt<-NA
 for (ag in 0:12)  w$wt[w$Age==ag]<-sw[,ag+1]
 
 w$wtp1<-0
 w$YC<-w$Year-w$Age
 for (i in 1:dim(w)[1]) try(w$wtp1[i]<-w$wt[w$YC==w$YC[i] & w$Age==w$Age[i]+1],silent=T)
 w$deltaW<-w$wtp1-w$wt
 
 # prepare the output object
 wt<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 wt[,ac(histPeriod),,,,] <-  t(sw[,-14])
 
 
 
 
 # fill in a 0 for age 0 
 wt[1,,,,,]<-0
 for( age in 1:12)
 {
 cat("age is",age,"\n")
 # model weight by an Arma process
 #first, age 1 
 dats<-w$wt[w$Age==age]
 try<-armaSearch(dats)
 coefs<-try$model
 p<-coefs[1]
 q<-coefs[2]
 ff<-paste("x~arma(",p,",",q,")",sep="")
 ff<-as.formula(ff)
 try<-armaFit( ff, data=dats)
 if (p==0){ aR<-c(0) } else {aR<-coef(try)[1:p]}
 if (q==0){ mA<-c(0) } else {mA<-coef(try)[(p+1):(p+q)]}
 
 for (i in 1:1000)
 {
 #rescale them
 x = armaSim(model = list(ar = aR, ma = mA), n = 1000)
 x<-x@.Data
 # standardise the series : same mean and var as the historical series
 x<-(x-mean(x))/sd(x)
 x<-0.8*sd(dats)*x+mean(dats)
 
 #remove the first simulated years until one is cloase enough to the last historical year to avoid a jump between historical and projection period
 last<-weighted.mean(dats[(length(dats)-2):length(dats)],(1:3)^2)
 xin<-(x>(last-0.015) & x<(0.015+last))
 xin<-cumsum(xin)
 x<-x[xin>0]
 x<-x[1:length(projPeriod)]
 x[x<min(dats)]<-min(dats)
 x<-c(dats,x)
 wt[age+1,,,,,i]<-x
 }
 }
swt<-stock.wt(stocks)
swt<-window(swt,start=1980,end=futureMaxYr)
swt<-propagate(swt,1000)
wt<-FLQuant(wt,dim=dim(swt),dimnames=dimnames(swt))

apply(wt,1:5,function(x) {sum(is.na(x))})

plot(wt)


save(wt,file=paste(inPath,"ARMAstock.wt.RData"))



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of catch weights 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#look at the difference between catch and stock weights
 sw<-read.table(file.path(inPath,"/",assess.name,"data","sw.dat"),skip=5)
 sw$Year<-1980:2013
 sw<-sw[sw$Year!=2013,]
 
 swt<-expand.grid(Year=1980:2012, Age=0:12)
 swt$swt<-NA
 for (ag in 0:12)  swt$swt[swt$Age==ag]<-sw[,ag+1]

 cw<-read.table(file.path(inPath,"/",assess.name,"data","cw.dat"),skip=5)
 cw$Year<-1980:2012
 
 w<-expand.grid(Year=1980:2012, Age=0:12)
 w$cwt<-NA
 for (ag in 0:12)  w$cwt[w$Age==ag]<-cw[,ag+1]

 w<-merge(swt,w)

 w$dev<-w$cwt-w$swt

library(lattice)
xyplot(dev~Year|Age,w)
histogram(~dev|Age,data=w)

devmu<-aggregate(dev~Age ,w,mean)
devsd<-aggregate(dev~Age ,w,sd)

# generate random deviations based on the historic ones, and add them to the stock weights to get the catch weights corresponding to the stock 
# weights. 
dev<-wt

  for(a in 0:12)  dev[ac(a),]<- rnorm(nits*dim(dev)[2],devmu$dev[devmu$Age==a],devsd$dev[devmu$Age==a])
cw<-wt+dev 
cw[,ac(histPeriod)]   <-  stocks@catch.wt[,ac(histPeriod)]

apply(cw,1:5,function(x) {sum(is.na(x))})
plot(cw)

save(cw,file=paste(inPath,"ARMAcatch.wt.RData"))
 
 
 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of maturity  option 1 : GLM model
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#load the data on maturity from WKPELA 2014
matPath<-"W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/stock weights/"

load(file=paste(matPath,"alldata.RData",sep=""))
dt<-dats[dats$origin=="catch" 
       & dats$age>0 
       & dats$age<12 
       & dats$month>=2 
       & dats$month<=7 
      # & dats$year>=1980 
       & !is.na(dats$age) 
       & !is.na(dats$month)
       ,]

dt$matur<-0
dt$matur[dt$mat!="i"]<-1

# fit a GLM model 
matnorm<-glm(matur~factor(age)+W,data=dt,family="binomial")


preddat<-expand.grid(age=1:11,W=50:1000)
preds<-predict(matnorm,newdata=preddat, type="response",se.fit=T)
preddat$pred<-preds[[1]]
preddat$predup<-preds[[1]]+1.96*preds[[2]]
preddat$predlo<-preds[[1]]-1.96*preds[[2]]

preddat<-preddat[preddat$age<9 & preddat$W<=400,]


# plot
image(1:8,50:400,(matrix(preddat$pred,nrow=8)), zlim=c(0,1), xlim=c(1,8), ylim=c(50,400),xlab="age",ylab="Weight",main="Probability of Being Mature (at a given age and weight)"  ,col = terrain.colors(12))
contour(1:8,50:400,(matrix(preddat$pred,nrow=8)), zlim=c(0,1), xlim=c(1,8), xlab="age",ylim=c(50,400),add=T,nlevels=3,levels=c(0.25,0.5,0.75),lwd=2)
contour(1:8,50:400,(matrix(preddat$predlo,nrow=8)), zlim=c(0,1), xlim=c(1,8), xlab="age",ylim=c(50,400),add=T,nlevels=3,levels=c(0.25,0.5,0.75),lty="dashed")
contour(1:8,50:400,(matrix(preddat$predup,nrow=8)), zlim=c(0,1), xlim=c(1,8), xlab="age",ylim=c(50,400),add=T,nlevels=3,levels=c(0.25,0.5,0.75),lty="dashed")

# predict maturity based on the catch weights
mat<-cw
mat["0",]<-0
mat[ac(4:12)] <- 1
for (a in 1:3) mat[ac(a),ac(projPeriod)] <-  apply(cw[ac(a),ac(projPeriod)],c(2:6),function(x) {predict(matnorm,newdata=data.frame(age=a,W=x*1000), type="response")})

mat[,ac(histPeriod)]   <-  stocks@mat[,ac(histPeriod)]

plot(mat)

save(mat,file=paste(inPath,"ARMAmatGLM.RData"))

                                                    # this doesn't really work this doesn't match with the historical data
                                                    

 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of maturity  option 1 : model as an ARMA process independent of the weights
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 mat<-read.table(file.path(inPath,"/",assess.name,"data","mo.dat"),skip=5)
 mat$Year<-1980:2013
 mat<-mat[mat$Year!=2013,]
 
 w<-expand.grid(Year=1980:2012, Age=0:12)
 w$mat<-NA
 for (ag in 0:12)  w$mat[w$Age==ag]<-mat[,ag+1]
 
 # prepare the output object
 mato<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 mato[,ac(histPeriod),,,,] <-  t(mat[,-14])
 
 # fill in a 0 for age 0 
 mato[1,,,,,]<-0
 mato[2,,,,,]<-0.11
 mato[5:13,,,,,]<-1
 for( age in 2:3)
 {
 cat("age is",age,"\n")
 # model weight by an Arma process
 #first, age 1 
 dats<-w$mat[w$Age==age]
 try<-armaSearch(dats)
 coefs<-try$model
 p<-coefs[1]
 q<-coefs[2]
 ff<-paste("x~arma(",p,",",q,")",sep="")
 ff<-as.formula(ff)
 try<-armaFit( ff, data=dats)
 if (p==0){ aR<-c(0) } else {aR<-coef(try)[1:p]}
 if (q==0){ mA<-c(0) } else {mA<-coef(try)[(p+1):(p+q)]}
 
 for (i in 1:1000)
 {
 #rescale them
 x = armaSim(model = list(ar = aR, ma = mA), n = 1000)
 x<-x@.Data
 # standardise the series : same mean and var as the historical series
 x<-(x-mean(x))/sd(x)
 x<-0.8*sd(dats)*x+mean(dats)
 
 #remove the first simulated years until one is cloase enough to the last historical year to avoid a jump between historical and projection period
 last<-weighted.mean(dats[(length(dats)-2):length(dats)],(1:3)^2)
 xin<-(x>(last-0.02) & x<(0.02+last))
 xin<-cumsum(xin)
 x<-x[xin>0]
 x<-x[1:length(projPeriod)]
 x[x<min(dats)]<-min(dats)
 x<-c(dats,x)
 mato[age+1,,,,,i]<-x
 }
 }
 
mat<-mat(stocks)
mat<-window(mat,start=1980,end=futureMaxYr)
mat<-propagate(mat,1000)
mat<-FLQuant(mato,dim=dim(mat),dimnames=dimnames(mat))

apply(mat,1:5,function(x) {sum(is.na(x))})
plot(mat)

save(mat,file=paste(inPath,"ARMAmat.RData"))



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of m prop
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 dt<-read.table(file.path(inPath,"/",assess.name,"data","pm.dat"),skip=5)
 dt$Year<-1980:2013
 dt<-dt[dt$Year!=2013,]
 
 w<-expand.grid(Year=1980:2012, Age=0:12)
 w$pm<-NA
 for (ag in 0:12)  w$pm[w$Age==ag]<-dt[,ag+1]
 
 # prepare the output object
 mp<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 mp[,ac(histPeriod),,,,] <-  t(dt[,-14])
 
 # fill in a 0 for age 0 
 
 age<-5
 dats<-w$pm[w$Age==age]
 try<-armaSearch(dats)
 coefs<-try$model
 p<-coefs[1]
 q<-coefs[2]
 ff<-paste("x~arma(",p,",",q,")",sep="")
 ff<-as.formula(ff)
 try<-armaFit( ff, data=dats)
 if (p==0){ aR<-c(0) } else {aR<-coef(try)[1:p]}
 if (q==0){ mA<-c(0) } else {mA<-coef(try)[(p+1):(p+q)]}
 
 for (i in 1:1000)
 {
 #rescale them
 x = armaSim(model = list(ar = aR, ma = mA), n = 10000)
 x<-x@.Data
 # standardise the series : same mean and var as the historical series
 x<-(x-mean(x))/sd(x)
 x<-0.8*sd(dats)*x+mean(dats)
 
 #remove the first simulated years until one is cloase enough to the last historical year to avoid a jump between historical and projection period
 last<-weighted.mean(dats[(length(dats)-2):length(dats)],(1:3)^2)
 xin<-(x>(last-0.01) & x<(0.01+last))
 xin<-cumsum(xin)
 x<-x[xin>0]
 x<-x[1:length(projPeriod)]
 x[x<min(dats)]<-min(dats)
 x<-c(dats,x)
mp[age+1,,,,,i]<-x
  }
 
for (i in 1:13)    mp[i,,,,,]   <-mp[age+1,,,,,]
 
 
mprop<-m.spwn(stocks)
mprop<-window(mprop,start=1980,end=futureMaxYr)
mprop<-propagate(mprop,1000)
mprop<-FLQuant(mp,dim=dim(mprop),dimnames=dimnames(mprop))

apply(mprop,1:5,function(x) {sum(is.na(x))})
plot(mprop)
save(mprop,file=paste(inPath,"ARMAmprop.RData"))


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   simulation of f prop using arma model  / computing it from the simulated m prop would require knowing (or simulating) how can are distributed in the year per age class.
#                                          / there is therefore no connection between changes in f and m prop here   
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


 dt<-read.table(file.path(inPath,"/",assess.name,"data","pf.dat"),skip=5)
 dt$Year<-1980:2013
 dt<-dt[dt$Year!=2013,]
 
 w<-expand.grid(Year=1980:2012, Age=0:12)
 w$pf<-NA
 for (ag in 0:12)  w$pf[w$Age==ag]<-dt[,ag+1]
 
 # prepare the output object
 fp<-array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,unit="unique",season="all",area="uniqu",iters=1:nits))
 fp[,ac(histPeriod),,,,] <-  t(dt[,-14])
 
 # fill in a 0 for age 0 
  fp[1,,,,,]<-0
 for( age in 1:12)
 {
 cat("age is",age,"\n")
 # model weight by an Arma process

 dats<-w$pf[w$Age==age]
 try<-armaSearch(dats)
 coefs<-try$model
 p<-coefs[1]
 q<-coefs[2]
 ff<-paste("x~arma(",p,",",q,")",sep="")
 ff<-as.formula(ff)
 try<-armaFit( ff, data=dats)
 if (p==0){ aR<-c(0) } else {aR<-coef(try)[1:p]}
 if (q==0){ mA<-c(0) } else {mA<-coef(try)[(p+1):(p+q)]}
 
 for (i in 1:1000)
 {
 #rescale them
 x = armaSim(model = list(ar = aR, ma = mA), n = 7000)
 x<-x@.Data
 # standardise the series : same mean and var as the historical series
 x<-(x-mean(x))/sd(x)
 x<-0.8*sd(dats)*x+mean(dats)
 
 #remove the first simulated years until one is cloase enough to the last historical year to avoid a jump between historical and projection period
 last<-weighted.mean(dats[(length(dats)-2):length(dats)],(1:3)^2)
 xin<-(x>(last-0.03) & x<(0.03+last))
 xin<-cumsum(xin)
 x<-x[xin>0]
 x<-x[1:length(projPeriod)]
 x[x<min(dats)]<-min(dats)
 x<-c(dats,x)
fp[age+1,,,,,i]<-x
 }
 }
 
 
fprop<-harvest.spwn(stocks)
fprop<-window(fprop,start=1980,end=futureMaxYr)
fprop<-propagate(fprop,1000)
fprop<-FLQuant(fp,dim=dim(fprop),dimnames=dimnames(fprop))

apply(fprop,1:5,function(x) {sum(is.na(x))})

plot(fprop)

save(fprop,file=paste(inPath,"ARMAfprop.RData"))



