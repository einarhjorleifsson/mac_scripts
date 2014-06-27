# RUN WITH R 2.13.2
rm(list=ls())


library(R2WinBUGS)
MyWinBugsDir <- "C:/Program Files (x86)/winbugs/WinBUGS14/"
library(FLCore)
library(FLAssess)


path          <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/"
inPath        <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/Data/"
codePath      <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/R code"
outPath       <- "W:/IMARES/Data/ICES-WG/WKPELA/WKPELA 2014 mackerel benchmark/MSE/exploration/plots/SR"
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/","/media/n/",path)
  inPath      <- sub("W:/","/media/n/",inPath)
  codePath    <- sub("W:/","/media/n/",codePath)
  outPath     <- sub("W:/","/media/n/",outPath)
}

  #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # 0): Read the data & assessment results (from the WKPELA 2014 benchmark assessment updated in April 2014 for the update advice)
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#- Read SAM output
setwd(path)

assess.name <-  "NEAMack-for-update-advice-2014"
source(paste(codePath,"/utils/flsam.r",sep=""))
run.dir<-paste(inPath,"/",assess.name,"/run",sep="")
source(paste(codePath,"/utils/SAM2FLR.r",sep=""))
name(Mac.sam) <- "NEA Mackerel"

recrPeriod<-1990:2011

# collect the SR pairs from the assessment output
Rtemp      <-  rec(Mac.sam)[,1:2]
SSBtemp    <-  ssb(Mac.sam)[,1:2]
names(Rtemp)<- c("year","R")
names(SSBtemp)<- c("year","SSB")
SRpairs<- merge(Rtemp,SSBtemp,all.y=T)
SRpairs<-SRpairs[is.element(SRpairs$year,as.numeric(recrPeriod)),]
rm(Rtemp,SSBtemp)


 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 1): Do the Bayesian estimation for the three models
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

set.seed(1) 
win.data<- list(n=dim(SRpairs)[1],
Recobs=log(SRpairs$R/100000),
SSB=SRpairs$SSB/1000000
)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Beverton and Holt model 
sink("WinBUGStest.txt")
cat("
model{
  #Diffuse Prior
BHa~dgamma(1.01,0.125)
BHb~dgamma(1.2,0.01)
tauBH~dgamma(0.001,0.001)
sigmaBH<-sqrt(1/tauBH)
 for(i in 1:n) {
       Recmod[i] <-log((SSB[i]/(BHa+SSB[i]*BHb)))
       Recobs[i] ~ dnorm(Recmod[i],tauBH)
}
}
",fill = TRUE)
sink()

#Set the initial values for the betas and sigma
inits <- function () {
   list(BHa=1,
BHb=5,
tauBH=.5)  }

#Parameters to estimate
params <- c("BHa","BHb","tauBH")


#MCMC settings
nc <- 3       #Number of chains
ni <- 1000     #Number of draws from posterior (for each chain)
nb <-  50     #Number of draws to discard as burn-in
nt <-   5     #Thinning rate

BevHolt <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "WinBUGStest.txt",
            n.thin = 5,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 40000,
            n.sims=100,
            debug = FALSE,
            bugs.directory = MyWinBugsDir)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Ricker model 
sink("WinBUGStest.txt")
cat("
model{
  #Diffuse Prior
BinfR~dgamma(1.01,0.125)
RecR~dgamma(1.5,.05)
tauR~dgamma(0.001,0.001)
sigmaR<-sqrt(1/tauR)

 for(i in 1:n) {
       Recmod[i] <-log(RecR *SSB[i]*exp(-BinfR*SSB[i]))
       Recobs[i] ~ dnorm(Recmod[i],tauR)   
}
}
",fill = TRUE)                
sink()

#Set the initial values for the betas and sigma
inits <- function () {
   list(BinfR=1e-4,
RecR=5,
tauR=.5)  }

#Parameters to estimate
params <- c("BinfR","RecR","tauR")


#MCMC settings

Ricker <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "WinBUGStest.txt",
            n.thin = 5,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 40000,
            n.sims=100,
            debug = F,
            bugs.directory = MyWinBugsDir)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Hockey stick model 
sink("WinBUGStest.txt")
cat("
model{
  #Diffuse Prior
HSb~dgamma(1.1,0.05)I(1.80465,4.16656)
#HSb~dgamma(1.1,1.25)I(1.80465,4.16656)
HSa~dgamma(1.3,0.05)
tauHS~dgamma(0.001,0.001)
sigmaHS<-sqrt(1/tauHS)
 for(i in 1:n) {
#          Recmod[i] <-log( HSa *step(SSB[i]-HSb)*HSb+HSa*SSB[i]*step(HSb-SSB[i]))
          Recmod[i] <-log( HSa *step(SSB[i]-HSb)+HSa*SSB[i]*step(HSb-SSB[i])/HSb)   # alternative formulation
        Recobs[i] ~ dnorm(Recmod[i],tauHS)
}
}
",fill = TRUE)
sink()

#Set the initial values for the betas and sigma
inits <- function () {
   list(HSb=2,
HSa=20,
tauHS=.5)  }

#Parameters to estimate
params <- c("HSb","HSa","tauHS")


#MCMC settings
nc <- 3       #Number of chains
ni <- 1000     #Number of draws from posterior (for each chain)
nb <-  50     #Number of draws to discard as burn-in
nt <-   5     #Thinning rate

HockeyStick <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "WinBUGStest.txt",
            n.thin = 5,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 40000,
            n.sims=100,
            debug = FALSE,
            bugs.directory = MyWinBugsDir)
            
            
            



 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 2): plot the models fitted ans the likelihood profiles
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


library(coda) # library just used for reading routines for winbugs
library(MASS)
pathFig<-paste(outPath,"/",sep="")

spp <- "MAC"
assMod <- "SAM"  # or "SCA".  Note: this needs to match the name given in the direcory where the data is.
srComb <-  "Three"#, "HS_RK" or "RK_BH"
resset<-"/"       # bayesian model fitted excluding 1982 YC
savePlots <- T           # save files or not
fileFormat <- "png"      # "eps", "png", "jpg"
if (!savePlots) dev.new()
options(graphics.record=T)
wth <- 7                                                                              
hght <- 7
dats<-SRpairs
 
### truncation values from separate analysis            
# (Half increment further out from max and min observed)
#       just tryied values a bit larger than the observed rct range / 100 000
Sru <- 150
Srl <- 5
  


#===========------------------- load data --------------------------===========#
#File names
Fname=c("BHa","BHb","BinfHS","BinfR","RecHS","RecR","sigmaBH","sigmaHS","sigmaR")

selection<-sample(c(1:length(BevHolt$sims.list$BHa)),2000,replace=F)
BHa<-BevHolt$sims.list$BHa[selection]
BHb<-BevHolt$sims.list$BHb[selection]
# Save orig.
BHasave=BHa; BHbsave=BHb
#Transform to the model spec used by FLR   
BHb=BHasave/BHbsave
BHa=1/BHbsave
sigmaBH<- (1/BevHolt$sims.list$tauBH[selection])^0.5


#   !!!! the HS formulation here is  R = A *SSB / B for SSB<B and R= A form SSB>B
HSa <-  HockeyStick$sims.list$HSa[selection] 
HSb <-  HockeyStick$sims.list$HSb[selection] 
sigmaHS<- (1/HockeyStick$sims.list$tauHS[selection])^0.5

#
Ra <-  Ricker$sims.list$RecR[selection] 
Rb <-  Ricker$sims.list$BinfR[selection] 
sigmaR<- (1/Ricker$sims.list$tauR[selection])^0.5



BHmods<-data.frame(BHa,BHb,sigmaBH)
Rikmods<-data.frame(Ra,Rb,sigmaR)
HSmods<-data.frame(HSa,HSb,sigmaHS)

###Stock recruit data       
# Must be same data used for Winbugs!!!
# Rec is the (Rec/ 1000), SSB is SSB/1000  .

Rec <- win.data$Recobs
SSB <- win.data$SSB
SSBO=SSB
RecO=Rec
Sssb=sort(SSB,index.return=TRUE)
SSB=Sssb$x
Rec=Rec[Sssb$ix]

#===========------------------ likelihoods -------------------------===========#
#### sections below calculate log likelihood of observations given set of models
### equations will depend on type of model
## here LL is LL fit to log of REC using Normal distribution formulated as 
#DM key part here - getting the probabilities for each model type (and using for proportions later on)

MCMCN=length(Ra)
SRpair=length(SSB)
Rsym=array(0,c(MCMCN,SRpair))
SSBp=seq(0,max(SSB),max(SSB)/50)

## Hockey Stick
PHS=0
P1=c(1:MCMCN)
for (i1 in 1:MCMCN) {
sll=0
sllq=0
for (i in 1:SRpair) {
  mu=log(HSa[i1]*(SSB[i]>=HSb[i1])+HSa[i1]*SSB[i]*(SSB[i]<HSb[i1])/HSb[i1]) # log Rec
  Rsym[i1,i]=mu
  ll=log(1/sqrt(2*3.14159*sigmaHS[i1]^2)*exp(-1/(2*(sigmaHS[i1]^2))*(mu-Rec[i])^2)) # LL of normal dist
  sll=sll+ll
}
P1[i1]=exp(sll)
}
PHS=sum(1/P1)/MCMCN
PHS=1/PHS

# Plot SR data and MCMC fits    
if (savePlots) png(file=paste(pathFig,spp,"_MCMC_SR-HS.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(SSB,exp(Rec), xlab="SSB", ylab="Recruitment",xlim=c(0,5),ylim=c(0,150))
for (i1 in seq(1,MCMCN,40)) lines(SSB,exp(Rsym[i1,]),col=5)
for (qnt in c(0.05,0.25,0.5,0.75,0.95)) lines(SSB,exp(apply(Rsym,2,quantile,qnt)),col=2)
lines(SSB,exp(Rsym[which.max(P1),]),col=1,lwd=2)
if (savePlots) dev.off()

## Ricker
PHR=0
P3=c(1:MCMCN)                                                                                                                                                          
for (i1 in 1:MCMCN) {                                                                                                                                            
sll=0                                                                                                                                                           
for (i in 1:SRpair) {                                                                                                                                               
  mu=log(Ra[i1]*SSB[i]*exp(-Rb[i1]*SSB[i]))
  Rsym[i1,i]=mu
  ll=log(1/sqrt(2*3.14159*sigmaR[i1]^2)*exp(-1/(2*(sigmaR[i1]^2))*(mu-Rec[i])^2))
  sll=sll+ll                                                                                                                                                    
}                                                                                                                                                               
P3[i1]=exp(sll)                                                                                                                                                      
}                                                                                                                                                               
PHR=sum(1/P3)/MCMCN                                                                                                                                              
PHR=1/PHR                                                                                                                                                       

# Plot SR data and MCMC fits    
if (savePlots) png(file=paste(pathFig,spp,"_MCMC_SR-Ric.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(SSB,exp(Rec), xlab="SSB", ylab="Recruitment",xlim=c(0,5),ylim=c(0,150))
for (i1 in seq(1,MCMCN,40)) lines(SSB,exp(Rsym[i1,]),col=5)
for (qnt in c(0.05,0.25,0.5,0.75,0.95)) lines(SSB,exp(apply(Rsym,2,quantile,qnt)),col=2)
lines(SSB,exp(Rsym[which.max(P3),]),col=1,lwd=2)
if (savePlots) dev.off()

##Bev Holt
PBH=0
P7=c(1:MCMCN)                                                                                                                                                          
for (i1 in 1:MCMCN) {                                                                                                                                            
sll=0                                                                                                                                                           
for (i in 1:SRpair) {
  mu <-log(BHa[i1]*SSB[i]/(BHb[i1]+SSB[i]))
  Rsym[i1,i]=mu
  ll=log(1/sqrt(2*3.14159*sigmaBH[i1]^2)*exp(-1/(2*(sigmaBH[i1]^2))*(mu-Rec[i])^2))
  sll=sll+ll                                                                                                                                                    

}                                                                                                                                                               
P7[i1]=exp(sll)                                                                                                                                                      
                                                                                                                          
}                                                                                                                                                               
PBH=sum(1/P7)/MCMCN                                                                                                                                              
PBH=1/PBH                                                                                                                                                       

# Plot SR data and MCMC fits  
if (savePlots) png(file=paste(pathFig,spp,"_MCMC_SR-BH.png",sep=""),height=hght,width=wth,units="in",res=200,)
plot(SSB,exp(Rec), xlab="SSB", ylab="Recruitment",xlim=c(0,5),ylim=c(0,150))
for (i1 in seq(1,MCMCN,40)) lines(SSB,exp(Rsym[i1,]),col=5)
for (qnt in c(0.05,0.25,0.5,0.75,0.95)) lines(SSB,exp(apply(Rsym,2,quantile,qnt)),col=2)
lines(SSB,exp(Rsym[which.max(P7),]),col=1,lwd=2)
if (savePlots) dev.off()




## prints out probabilities of each of the three models picked up   
pp=seq(1,MCMCN,1)
P1a=sort(P1)[pp]
P3a=sort(P3)[pp]
P7a=sort(P7)[pp]

Total=1
for (nn in seq(1,500,2)) {
pp=seq(nn,MCMCN,1)
P1a=sort(P1/Total)[pp]
P3a=sort(P3/Total)[pp]
P7a=sort(P7/Total)[pp]
P1t=1/mean(1/P1a)
P3t=1/mean(1/P3a)
P7t=1/mean(1/P7a)

print(c(P1t,P3t)/sum(P1t,P3t))
}

pp=seq(1,MCMCN,1)
P1a=sort(P1/Total)[pp]
P3a=sort(P3/Total)[pp]
P7a=sort(P7/Total)[pp]

# Plot model likelihoods
if (savePlots) png(file=paste(pathFig,spp,"_Bayes-Lkhd.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(x=pp,y=log(P1a),type="n",main="Likelihood of Bayes model",xlab="model order",ylab="log likelihood")#,ylim=c(-50,-36)) #SOL ylim=c(-75,-55)
lines(x=pp,y=log(P1a),lty=1,col=1)
lines(x=pp,y=log(P3a),lty=2,col=2)
lines(x=pp,y=log(P7a),lty=3,col=3)
legend(x="bottomright",legend=c("HS","Rk","BH"),lty=c(1,2,3),col=c(1,2,3))
if (savePlots) dev.off()

##DM getting coefficients for max likelihood model fit
pp1=which.max(P1)             
pp3=which.max(P3)
pp7=which.max(P7)
maxLFit <- paste('H-Stick : A=',HSa[pp1]," B=",HSb[pp1]," sigma=",sigmaHS[pp1],
';Ricker  : A=',Ra[pp3]," B=",Rb[pp3]," sigma=",sigmaR[pp3],
';Bev-Holt: A=',BHa[pp7]," B=",BHb[pp7]," sigma=",sigmaBH[pp7], sep=" ") 

write(maxLFit, file=paste(pathFig,spp,"_MaxLikelihoodFits.txt",sep=""),sep=",")



#  plot fitted vs. observed recruitment
 png(file=paste(pathFig,spp,"fitted vs observed.png",sep=""),height=6,width=12,units="in",res=200)
 par(mfrow=c(1,3))
# for Ricker
mu=(Ra[pp3]*SSB*exp(-Rb[pp3]*SSB))
plot(exp(Rec),mu, xlim=c(0,150),ylim=c(0,75),xlab="observer Rct",ylab="fitted",main="Ricker")
abline(0,1)
mu=(HSa[pp1]*(SSB>=HSb[pp1])+HSa[pp1]*SSB*(SSB<HSb[pp1])/HSb[pp1])
plot(exp(Rec),mu, xlim=c(0,150),ylim=c(0,75),xlab="observer Rct",ylab="fitted",main="Hockey Stick")
abline(0,1)
mu <-(BHa[pp7]*SSB/(BHb[pp7]+SSB))
plot(exp(Rec),mu, xlim=c(0,150),ylim=c(0,75),xlab="observer Rct",ylab="fitted",main="Beverton and Holt")
abline(0,1)
dev.off()




 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 3): probabilities for each model
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# change start values to get correct total proportions
srComb<-"Three"
if (srComb=="Three") {
  # three models
  srProbs <- c(PHS,PHR,PBH)/sum(PHS,PHR,PBH)
  probHS <- srProbs[1]
  probHR <- srProbs[2]
  probBH <- srProbs[3]
  # HS RK BH = 0.1421944 0.8578056 0.0000000         #DM these values are hard wired
  mod=array("HSL",1000)
  for (i in round(seq(2,1000,(probHS+probHR)/probHR))){mod[i]="RKL"}
  for (i in round(seq(3.4,1000,1/probBH))){mod[i]="BHL"}
  #DM check to see if ended up with right numbers
  print(c(sum(mod=="HSL"),sum(mod=="RKL"),sum(mod=="BHL")))           
  
  } else if (srComb=="HS_BH") {
  # two models: HS and Ricker
  srProbs <- c(PHS,PBH)/sum(PHS,PBH)
  probHS <- srProbs[1]
  probBH <- srProbs[2]
  # HS RK = 0.5832516	0.4167484
  mod=array("HSL",1000)
  for (i in round(seq(1.5,1000,1/probBH))){mod[i]="BHL"}
  #DM check to see if ended up with right numbers
  print(c(sum(mod=="HSL"),sum(mod=="BHL")))           
  
  } else if (srComb=="RK_BH") {
  # two models: Ricker and BevHolt
  srProbs <- c(PHR,PBH)/sum(PHR,PBH)
  probHR <- srProbs[1]
  probBH <- srProbs[2]
  # HS RK = 0.5832516	0.4167484
  mod=array("RKL",1000)
  for (i in round(seq(1.5,1000,1/probBH))){mod[i]="BHL"}
  #DM check to see if ended up with right numbers
  print(c(sum(mod=="RKL"),sum(mod=="BHL")))           
  }


 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 4): otholith plot of the parameter estimates for the three models
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

selection<-c(length(BevHolt$sims.list$BHa):1)[1:20000]

# Loop over SR types
for (sr in c("HS","RH","BH")) {

#First save max LL parameter values from input data
if (sr=="HS") {
# log lilihood vbalues for HS 
LLa = HSa[pp1]
LLb = HSb[pp1]
LLsigma=sigmaHS[pp1]
HSa <-  HockeyStick$sims.list$HSa[selection] 
HSb <-  HockeyStick$sims.list$HSb[selection] 
sigmaHS<- (1/HockeyStick$sims.list$tauHS[selection])^0.5
AA=HSa
BB=HSb

} else if (sr=="RH") {
# log lilihood vbalues for Ricker 
LLa = Ra[pp3]
LLb = Rb[pp3]
LLsigma=sigmaR[pp3]
Ra <-  Ricker$sims.list$RecR[selection] 
Rb <-  Ricker$sims.list$BinfR[selection] 
sigmaR<- (1/Ricker$sims.list$tauR[selection])^0.5
AA=Ra
BB=Rb

} else if (sr=="BH") {

LLa = BHasave[pp7]
LLb = BHbsave[pp7]
LLsigma=sigmaBH[pp7]
BHa<-BevHolt$sims.list$BHa[selection]
BHb<-BevHolt$sims.list$BHb[selection]
# Save orig.
BHasave=BHa; BHbsave=BHb
BHb=BHasave/BHbsave
BHa=1/BHbsave
sigmaBH<- (1/BevHolt$sims.list$tauBH[selection])^0.5
AA=BHasave
BB=BHbsave
}

### run this plotiing bit for each bayes model ------------ Routine borrowed and amended from Mark Payne
AA.est=median(AA)
BB.est=median(BB)

# some settings normally set in the subroutine call
n.grid=50    # this and the n below will control smoothness
show.points=TRUE
show.ll=TRUE
do.contours=TRUE
filled.contours=FALSE
f.ages=NULL
margin.plots=TRUE
xlim= quantile(BB,0.99)
ylim= quantile(AA,0.99)
debug=FALSE
n=20000 # set to match my full data length
pch="."
show.grid=TRUE
alpha=0.05   # sets intervals on pfs
show.estimate=TRUE
thin=20
TH=seq(thin/2,n,thin)

# main contouring section

    kern  <-  kde2d(BB,AA,n=n.grid)
    #Calculate cumulative distribution function
    kz    <-  as.vector(kern$z)
    ord   <-  order(kz)
    cumfrac <-  cumsum(kz[ord])/sum(kz)
    cumfrac.matrix  <-  matrix(cumfrac[rank(kz)],nrow=nrow(kern$z),ncol=ncol(kern$z))
    contour.args <- list()
    contour.args$levels <-   c(0.1,0.25,0.50,0.75,0.90)
    contour.args$labels <-   NULL
#      if(is.null(contour.args$lty))   contour.args$lty <-   c(1,1,2,3,4)
#      if(is.null(contour.args$lwd))   contour.args$lwd <-   c(1,3,1,1,1)
    if(is.null(contour.args$method))contour.args$method <-    "edge"
#    if(is.null(contour.args$labels))contour.args$labels <-    NULL
#      if(filled.contours) {      
#       do.call(filled.contour,c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,nlevels=100,color.palette=heat.colors))
#    }
    otolith.obj  <-  c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,contour.args)
# this form works but gives odd contour labeling  - alternative crude seting at end to deal with issue
# I cannot work out how to change contour.args the iff line 6 above commented out was an attemp that did not work 

if (savePlots) png(file=paste(pathFig,spp,"_",sr,"_Otolith.png",sep=""),height=hght,width=wth,units="in",res=200)

    if(!show.points) { pch <- NA }
    if(margin.plots) {
      layout(matrix(c(1,4,3,2),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
      par(mar=c(0.5,0.5,0.5,0.5),oma=c(5.1,4.1,4.1,2.1))
    }


#    x.lab <-  "Q" # not used I think
   xlim  <- range(pretty(BB[BB<quantile(BB,0.99)]))
   ylim  <- range(pretty(AA[AA<quantile(AA,0.99)]))
    #First the horizontal plot
    if(margin.plots) { 
      densF   <-  density(BB)
      plot(densF,ann=FALSE,xaxt="n",yaxt="n",type="l",xlim=xlim)   
      if(show.grid) grid()      
      title(ylab="Prob. Density",xpd=NA,mgp=c(1,1,0))
      #Calculate 95% confidence intervals
      cumsumF.fun <-  approxfun(cumsum(densF$y)/sum(densF$y),densF$x)
      densF.fun   <-  approxfun(densF$x,densF$y)
      ul.F    <-  cumsumF.fun(1-alpha/2)      
      ul.dens <-  densF.fun(ul.F)
      ll.F    <-  cumsumF.fun(alpha/2)      
      ll.dens <-  densF.fun(ll.F)
      points(c(ll.F,ul.F),c(ll.dens,ul.dens),pch="|",cex=1.5)
      text(c(ll.F,ul.F),c(ll.dens,ul.dens),label=sprintf("%.3f",round(c(ll.F,ul.F),3)),pos=4)
      if(show.estimate) { 
        points(BB.est,densF.fun(BB.est),pch=19,cex=1.5)
        text(BB.est,densF.fun(BB.est),label=sprintf("%.3f",round(BB.est,3)),pos=4)
        }
      if(show.ll) { 
        points(LLb,densF.fun(LLb),pch=19,cex=1.5,col=4)
        text(LLb,densF.fun(LLb),label=sprintf("%.3f",round(LLb,3)),pos=4)
      }
    }
    #Now the vertical plot
    if(margin.plots) { 
      densAA <-  density(AA)
      plot(densAA$y,densAA$x,xaxt="n",yaxt="n",type="l",ylim=ylim)
      abline(v=0,col="grey")
      if(show.grid) grid()
      title(xlab="Prob. Density",xpd=NA,mgp=c(1,1,0))      
      #Calculate 95% confidence intervals
      cumsumAA.fun <-  approxfun(cumsum(densAA$y)/sum(densAA$y),densAA$x)
      densAA.fun   <-  approxfun(densAA$x,densAA$y)
      ul.AA    <-  cumsumAA.fun(1-alpha/2)      
      ul.dens <-  densAA.fun(ul.AA)
      ll.AA    <-  cumsumAA.fun(alpha/2)      
      ll.dens <-  densAA.fun(ll.AA)
      points(c(ll.dens,ul.dens),c(ll.AA,ul.AA),pch="-",cex=2)
      text(c(ll.dens,ul.dens),c(ll.AA,ul.AA),label=round(c(ll.AA,ul.AA),3),pos=4)
      if(show.estimate) {
        points(densAA.fun(AA.est),AA.est,pch=19,cex=1.5)      
        text(densAA.fun(AA.est),AA.est,,label=round(AA.est,3),pos=2)
        }
      if(show.ll) { 
        points(densAA.fun(LLa),LLa,pch=19,cex=1.5,col=4)      
        text(densAA.fun(LLa),LLa,,label=round(LLa,3),pos=2)
      }
    }
    #Now the main plot
    plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="")
#   Three sets of points for 3 Bayes MCMC chains - normally one would do
    if(show.points) points(BB[TH],AA[TH],pch=19,col=2,cex=0.05)
    title(xlab="B",ylab="A",xpd=NA)
    if(show.estimate) points(BB.est,AA.est,pch=19,cex=1.5)
    if(show.ll) points(LLb,LLa,pch=19,col=4,cex=1.5)
    if(show.grid) grid()
    contour(otolith.obj,levels=contour.args$levels,lables=NULL,add=TRUE)
# dont run lines below
    if(do.contours) {
      do.call(contour,otolith.obj)
    }
if (savePlots) dev.off()
  } 
  


 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 5): computing the percentage of points outside the limits for authorised recruitment values
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  

P1=c(1:2000)
Sr=array(0,c(SRpair,3000))
dev=array(0,c(SRpair,3000))
HSNmu=array(0,c(SRpair,2000))
lowlim=array(0,2000)
uplim=array(0,2000)
RecHSNub=array(0,2000)
mincnt=array(0,2000)
maxcnt=array(0,2000)
RecHSNmn=array(0,2000)
for (i1 in 1:2000) {
sll=0
for (i in 1:SRpair) {
  mu=log(HSa[i1]*(SSB[i]>=HSb[i1])+HSa[i1]*SSB[i]*(SSB[i]<HSb[i1])/HSb[i1])
  dev[i,]=rnorm(3000,0,sigmaHS[i1])
  Sr[i,]=exp(mu+dev[i,])
  HSNmu[i,i1]=mean(Sr[i,])
}
mincnt[i1]=sum(Sr<Srl)
maxcnt[i1]=sum(Sr>Sru)
lowlim[i1] = min(dev[(Sr>Srl)])
uplim[i1] = max(dev[(Sr<Sru)]) 


RecHSNub[i1]=HSa[i1]+log(mean(HSNmu[,i1])/mean(Sr*(Sr>Srl)*(Sr<Sru)))
RecHSNmn[i1]=mean(HSNmu[,i1])
P1[i1]=exp(sll)
}
PHS=sum(1/P1)/2000
PHS=1/PHS
mincnt=mincnt/SRpair/3000
maxcnt=maxcnt/SRpair/3000
HSN=as.data.frame(cbind(HSa,HSb,sigmaHS,mincnt,maxcnt,lowlim,uplim,RecHSNmn,RecHSNub))



PHR=0
P2=c(1:2000)                                                                                                                                                          
RNmu=array(0,c(SRpair,2000))
RecRNub=array(0,2000)
RecRNmn=array(0,2000)
for (i1 in 1:2000) {
sll=0                                                                                                                                                           
for (i in 1:SRpair) {
  mu=log(Ra[i1]*SSB[i]*exp(-Rb[i1]*SSB[i]))
  dev[i,]=rnorm(3000,0,sigmaR[i1])
  Sr[i,]=exp(mu+dev[i,])
  RNmu[i,i1]=mean(Sr[i,])
}
mincnt[i1]=sum(Sr<Srl)
maxcnt[i1]=sum(Sr>Sru)
lowlim[i1] = min(dev[(Sr>Srl)])
uplim[i1] = max(dev[(Sr<Sru)]) 

RecRNub[i1]=Ra[i1]+log(mean(RNmu[,i1])/mean(Sr*(Sr>Srl)*(Sr<Sru)))
RecRNmn[i1]=mean(RNmu[,i1])

P2[i1]=exp(sll)
}                                                                                                                                                               
PHR=sum(1/P2)/2000                                                                                                                                              
PHR=1/PHR                                                                                                                                                       
mincnt=mincnt/3000/SRpair
maxcnt=maxcnt/3000/SRpair
RkN=as.data.frame(cbind(Ra,Rb,sigmaR,mincnt,maxcnt,lowlim,uplim,RecRNmn,RecRNub))


BHN=0
P4=c(1:2000)                                                                                                                                                          
BHNmu=array(0,c(SRpair,2000))
RecBHNub=array(0,2000)
RecBHNmn=array(0,2000)
for (i1 in 1:2000) {
for (i in 1:SRpair) {
  mu <-log(BHa[i1]*SSB[i]/(BHb[i1]+SSB[i]))
  dev[i,]=rnorm(3000,0,sigmaBH[i1])
  Sr[i,]=exp(mu+dev[i,])
  BHNmu[i,i1]=mean(Sr[i,])
#print(c(x,log(Rec[i])))
}                                                                                                                                                               
mincnt[i1]=sum(Sr<Srl)
maxcnt[i1]=sum(Sr>Sru)
lowlim[i1] = min(dev[(Sr>Srl)])
uplim[i1] = max(dev[(Sr<Sru)]) 

RecBHNub[i1]=BHa[i1]+log(mean(BHNmu[,i1])/mean(Sr*(Sr>Srl)*(Sr<Sru)))
RecBHNmn[i1]=mean(BHNmu[,i1])


}                                                                                                                                                               
mincnt=mincnt/3000/SRpair
maxcnt=maxcnt/3000/SRpair
BHN=as.data.frame(cbind(BHa,BHb,sigmaBH,mincnt,maxcnt,lowlim,uplim,RecBHNmn,RecBHNub))

print(c(sum(mod=="RKL"),sum(mod=="HSL"),sum(mod=="BHL")))







 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 6): prepare models sets with and without truncation
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# without tuncation

pp=c(1:1000)
HSLp=pp[mod=="HSL"]
RKLp=pp[mod=="RKL"]
BHLp=pp[mod=="BHL"]

untrncModset=as.data.frame(cbind(mod,HSN[1:1000,1],HSN[1:1000,2:3],HSN[1:1000,6:7],HSN[1:1000,4:5],HSN[1:1000,9]/HSN[1:1000,1]))
untrncModset[,-1]<-NA
i=RKLp
untrncModset[i,2]=RkN[i,1]
untrncModset[i,3:4]=RkN[i,2:3]
untrncModset[i,5:6]=RkN[i,6:7]
untrncModset[i,7:8]=RkN[i,4:5]
untrncModset[i,9]=RkN[i,9]/RkN[i,1]

i=BHLp
untrncModset[i,2]=BHN[i,1]
untrncModset[i,3:4]=BHN[i,2:3]
untrncModset[i,5:6]=BHN[i,6:7]
untrncModset[i,7:8]=BHN[i,4:5]
untrncModset[i,9]=BHN[i,9]/BHN[i,1]
#
i=HSLp
untrncModset[i,2]=HSN[i,1]
untrncModset[i,3:4]=HSN[i,2:3]
untrncModset[i,5:6]=HSN[i,6:7]
untrncModset[i,7:8]=HSN[i,4:5]
untrncModset[i,9]=HSN[i,9]/HSN[i,1]

names(untrncModset)[2:4]=c("A","B","sigma")
names(untrncModset)[9]="Cfac"

if (savePlots) png(file=paste(pathFig,spp,"_ModsUpperlim_unTrunc.png",sep=""),height=hght,width=wth,units="in",res=200)
hist(untrncModset[,8]*100,30,main="No of Models exceeding upper limit", xlab="Percentage of simulated values exceeding upper limit" )
if (savePlots) dev.off()

if (savePlots) png(file=paste(pathFig,spp,"_ModsLowerlim_unTrunc.png",sep=""),height=hght,width=wth,units="in",res=200)
hist(untrncModset[,7]*100,30,main="No of Models exceeding lower limit", xlab="Percentage of simulated values exceeding lower limit" )
if (savePlots) dev.off()



# with tuncation

pp=c(1:1000)
HSLp=pp[mod=="HSL"]
RKLp=pp[mod=="RKL"]
BHLp=pp[mod=="BHL"]

trncModset=as.data.frame(cbind(mod,HSN[1:1000,9],HSN[1:1000,2:3],HSN[1:1000,6:7],HSN[1:1000,4:5],HSN[1:1000,9]/HSN[1:1000,1]))

i=RKLp
trncModset[i,2]=RkN[i,9]
trncModset[i,3:4]=RkN[i,2:3]
trncModset[i,5:6]=RkN[i,6:7]
trncModset[i,7:8]=RkN[i,4:5]
trncModset[i,9]=RkN[i,9]/RkN[i,1]

i=BHLp
trncModset[i,2]=BHN[i,9]
trncModset[i,3:4]=BHN[i,2:3]
trncModset[i,5:6]=BHN[i,6:7]
trncModset[i,7:8]=BHN[i,4:5]
trncModset[i,9]=BHN[i,9]/BHN[i,1]
#
i=HSLp
trncModset[i,2]=HSN[i,9]
trncModset[i,3:4]=HSN[i,2:3]
trncModset[i,5:6]=HSN[i,6:7]
trncModset[i,7:8]=HSN[i,4:5]
trncModset[i,9]=HSN[i,9]/HSN[i,1]


names(trncModset)[2:4]=c("A","B","sigma")
names(trncModset)[9]="Cfac"

if (savePlots) png(file=paste(pathFig,spp,"_ModsUpperlim_trunc.png",sep=""),height=hght,width=wth,units="in",res=200)
hist(trncModset[,8]*100,30,main="No of Models exceeding upper limit", xlab="Percentage of simulated values exceeding upper limit" )
if (savePlots) dev.off()

if (savePlots) png(file=paste(pathFig,spp,"_ModsLowerlim_trunc.png",sep=""),height=hght,width=wth,units="in",res=200)
hist(trncModset[,7]*100,30,main="No of Models exceeding lower limit", xlab="Percentage of simulated values exceeding lower limit" )
if (savePlots) dev.off()

 #end of not truncated loop
      
#########################################################
#DM saving the model parameters
write.table(trncModset[,1:6],file=paste(pathFig,spp,"_modelparams_trnc.dat",sep=""),row.names=FALSE,col.names=TRUE,sep=" ")
write.table(untrncModset[,1:6],file=paste(pathFig,spp,"_modelparams.dat",sep=""),row.names=FALSE,col.names=TRUE,sep=" ")
write.table(trncModset,file=paste(pathFig,spp,"_modelparamsfull_trnc.dat",sep="_"),row.names=FALSE,col.names=TRUE,sep=" ")
write.table(untrncModset,file=paste(pathFig,spp,"_modelparamsfull.dat",sep="_"),row.names=FALSE,col.names=TRUE,sep=" ")




 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 7): plot simulated recruitment based on the truncated and untruncated model sets
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

trnc<-F
if(trnc) Modset <- trncModset else Modset <- untrncModset

# plot out recruits
dats<-SRpairs
Rec <- log(dats$R/100000)
SSB <- dats$SSB/1e6

SSBO=SSB
RecO=Rec
Sssb=sort(SSB,index.return=TRUE)
SSB=Sssb$x
Rec=Rec[Sssb$ix]
mn=length(RecO)
SSBs=seq(1:mn*40*1000)
Recs=seq(1:mn*40*1000)

### Untruncated

if (savePlots) png(file=paste(pathFig,"_MergedSRfit_unTrunc.png",sep=""),height=hght,width=wth,units="in",res=200)
repl<-1
Xmax <- ceiling(max(SSBO, na.rm=T))
Ymax <- ceiling(max(exp(RecO), na.rm=T))
plot(SSBO,exp(RecO),xlim=c(0,Xmax),ylim=c(0,150),type="p",pch=19,col=10,xlab="SSB (millions of tonnes)",ylab="Recruits (*1E+8)", main=spp)
abline(h=0)
for (i in 1:1000) {
#SSB = runif(3000,0.184,1.171)
#SSB = seq(0,max(SSBO),length.out=length(SSBO)*repl) ## to match SSB values to original
SSB = rep(SSBO,repl) ## to match SSB values to original
 
 if (Modset[i,1]=="HSL") {
  mu=log(Modset$A[i]*(SSB>=Modset$B[i])+Modset$A[i]*SSB*(SSB<Modset$B[i])/Modset$B[i])
  R2=rnorm(mn*repl,0,Modset$sigma[i])
  Rec=exp(mu+R2)
 }
 if (Modset[i,1]=="RKL") {
  mu=log(Modset$A[i]*SSB*exp(-Modset$B[i]*SSB))
  R2=rnorm(mn*repl,0,Modset$sigma[i])
  Rec=exp(mu+R2[1:mn*repl])
 }
 if (Modset[i,1]=="BHL") {
  mu <-log(Modset$A[i]*SSB/(Modset$B[i]+SSB))
  R2=rnorm(mn*repl,0,Modset$sigma[i])
  Rec=exp(mu+R2[1:mn*repl])
 }
points(SSB,Rec,type="p",pch=20,col=1,cex=0.0625)
SSBs[((i-1)*mn*repl+1):(i*mn*repl)]=SSB
Recs[((i-1)*mn*repl+1):(i*mn*repl)]=Rec
}
Rsims<-Recs
Correct<-data.frame(trnc=c(F,T),Correct=c(1,1))
Correct[1,2]=mean(Recs)/mean(exp(RecO))
points(SSBO,exp(RecO),type="p",pch=19,col=10,cex=1.25)


seqMin <- 0.5
seqStep <- 0.5
Xmax <-5.5
up=seq(seqMin,Xmax,seqStep)
lw=seq(seqMin,Xmax,seqStep)
md=seq(seqMin,Xmax,seqStep)
ssb=seq(seqMin,Xmax,seqStep)
loopNum <- length(ssb)    
for (j in 1:loopNum){     
  up[j]=quantile(Recs[((SSBs>up[j])*(SSBs<up[j+1]))>0],probs=.95,na.rm=TRUE)
  lw[j]=quantile(Recs[((SSBs>lw[j])&(SSBs<lw[j+1]))>0],probs=.05,na.rm=TRUE)
  md[j]=quantile(Recs[((SSBs>md[j])&(SSBs<md[j+1]))>0],probs=.5,na.rm=TRUE)
  }
  
lines(ssb,md[1:loopNum],col="green",lwd=3)
lines(ssb,up[1:loopNum],col=4,lwd=3)
lines(ssb,lw[1:loopNum],col=4,lwd=3)

if (savePlots) dev.off()


# plot out recruits  with truncation
trnc<-T
if(trnc) Modset <- trncModset else Modset <- untrncModset
if (savePlots) png(file=paste(pathFig,spp,"_MergedSRfit_trunc.png",sep=""),height=hght,width=wth,units="in",res=200)
Xmax <- ceiling(max(SSBO, na.rm=T))
Ymax <- ceiling(max(exp(RecO), na.rm=T))
plot(SSBO,exp(RecO),xlim=c(0,Xmax),ylim=c(0,150),type="p",pch=19,col=10,xlab="SSB (millions of tonnes)",ylab="Recruits (*1E8)", main=spp)
for (i in 1:1000) {
#SSB = runif(mn*40,.1,5.0)
SSB = rep(SSBO,40) ## to match SSB values to original
 if (Modset[i,1]=="HSL") {
  mu=log(Modset$A[i]*(SSB>=Modset$B[i])+Modset$A[i]*SSB*(SSB<Modset$B[i])/Modset$B[i])
  R1=rnorm(4000,0,Modset$sigma[i])
  R2=R1[R1>Modset$lowlim[i]]
  R2=R2[R2<Modset$uplim[i]]
  Rec=exp(mu+R2[1:mn*40])
 }
 if (Modset[i,1]=="RKL") {
  mu=log(Modset$A[i]*SSB*exp(-Modset$B[i]*SSB))
  R1=rnorm(4080,0,Modset$sigma[i])  #SOL R1=rnorm(4000,0,Modset$sigma[i])
  R2=R1[R1>Modset$lowlim[i]]
  R2=R2[R2<Modset$uplim[i]]
  Rec=exp(mu+R2[1:mn*40])
 }
 if (Modset[i,1]=="BHL") {
  mu <-log(Modset$A[i]*SSB/(Modset$B[i]+SSB))
  R1=rnorm(4080,0,Modset$sigma[i]) #SOL R1=rnorm(4000,0,Modset$sigma[i])
  R2=R1[R1>Modset$lowlim[i]]
  R2=R2[R2<Modset$uplim[i]]
  Rec=exp(mu+R2[1:mn*40])
 }
points(SSB,Rec,type="p",pch=20,col=1,cex=0.0625)
SSBs[((i-1)*mn*40+1):(i*mn*40)]=SSB
Recs[((i-1)*mn*40+1):(i*mn*40)]=Rec
}
points(SSBO,exp(RecO),type="p",pch=19,col=10,cex=1.25)


#PLE
if (spp=="HOM") {
  seqMin <- 0.1
  seqStep <- 0.5
#  loopNum <- 24    
  }

up=seq(seqMin,Xmax,seqStep)
lw=seq(seqMin,Xmax,seqStep)
md=seq(seqMin,Xmax,seqStep)
ssb=seq(seqMin,Xmax,seqStep)
loopNum <- length(ssb)    
for (j in 1:loopNum){     
  up[j]=quantile(Recs[((SSBs>up[j])*(SSBs<up[j+1]))>0],probs=.95,na.rm=TRUE)
  lw[j]=quantile(Recs[((SSBs>lw[j])&(SSBs<lw[j+1]))>0],probs=.05,na.rm=TRUE)
  md[j]=quantile(Recs[((SSBs>md[j])&(SSBs<md[j+1]))>0],probs=.5,na.rm=TRUE)
  }
  
lines(ssb,md[1:loopNum],col=7,lwd=3)
lines(ssb,up[1:loopNum],col=4,lwd=3)
lines(ssb,lw[1:loopNum],col=4,lwd=3)

if (savePlots) dev.off()
#DM ?

Rsims<-Recs
Correct[2,2]=mean(Recs)/mean(exp(RecO))
#DM writing out corrected values
write.table(Modset[,1:6],file=paste(spp,"_modelparamsCor",if(trnc) "_trnc",".dat",sep=""),row.names=FALSE,col.names=TRUE,sep=" ")
write.table(Modset,file=paste(spp,"_modelparamsCorfull",if(trnc) "_trnc",".dat",sep="_"),row.names=FALSE,col.names=TRUE,sep=" ")





 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # 7): select the final model set and rescale the parameters to the correct dimention
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# final choice about truncation
trnc<-F
if(trnc) Modset <- trncModset else Modset <- untrncModset


# apply a correction to make sure the simulated mean is the same as the data mean
Core<-Correct$Correct[Correct$trnc==trnc]
Modset$A=Modset$A/Correct 


# rescale the parameters
Modset<-Modset[,1:4]
Modset$A[Modset$mod=="HSL"]<-Modset$A[Modset$mod=="HSL"]*1e5
Modset$B[Modset$mod=="HSL"]<-Modset$B[Modset$mod=="HSL"]*1e6

Modset$A[Modset$mod=="BHL"]<-Modset$A[Modset$mod=="BHL"]*1e5
Modset$B[Modset$mod=="BHL"]<-Modset$B[Modset$mod=="BHL"]*1e6

Modset$A[Modset$mod=="RKL"]<-Modset$A[Modset$mod=="RKL"]/10
Modset$B[Modset$mod=="RKL"]<-Modset$B[Modset$mod=="RKL"]*1e-6


#final check that there is no mistake
a<-aggregate(Modset$A,by=list(Modset$mod),FUN=median)
b<-aggregate(Modset$B,by=list(Modset$mod),FUN=median)

ssB<-seq(0,max(dats$SSB),length.out=50)

plot(dats$SSB,dats$R,xlim=c(0,5e6),ylim=c(0,1.1e7))

rk<-a[3,2]*ssB*exp(-ssB*b[3,2])
lines(ssB,rk)

hs<-(a[2,2]*ssB*(ssB<b[2,2])/b[2,2] ) + a[2,2] *(ssB>=b[2,2])
lines(ssB,hs)

bh<-a[1,2]*ssB/(b[1,2]+ssB)
lines(ssB,bh)


Modset$mod[Modset$mod=="HSL"]<-"segreg"
Modset$mod[Modset$mod=="BHL"]<-"bevholt"
Modset$mod[Modset$mod=="RKL"]<-"ricker"

write.csv(Modset,file=paste(inPath,"SRbayes",min(recrPeriod),"_",max(recrPeriod),".csv",sep=""))

