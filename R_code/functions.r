#-------------------------------------------------------------------------------
#- Conversions
#-------------------------------------------------------------------------------
an <- function(x){return(as.numeric(x))}
ac <- function(x){return(as.character(x))}



#-------------------------------------------------------------------------------
#- simple plot function for checking input
#-------------------------------------------------------------------------------

plquant<-function(quant)
{
library(lattice)
xyplot(data~year|age,groups=iter,data=quant,type="l")
}

#-------------------------------------------------------------------------------
#- Draw random numbers from a log-normal distribution, but with truncation possibility
#-------------------------------------------------------------------------------
rtlnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
{
   ret <- numeric()
   if (length(n) > 1)
       n <- length(n)
   while (length(ret) < n) {
       y <- rlnorm(n - length(ret), mean, sd)
       y <- y[y >= lower & y <= upper]
       ret <- c(ret, y)
   }
   stopifnot(length(ret) == n)
   ret
}

#-------------------------------------------------------------------------------
#- Fit vonBertalanffy curve to weight data and draw new values based on sd
#-------------------------------------------------------------------------------
randomWeightBertalanffy <- function(histWeight,ages,histPeriod,projPeriod,iters,estim=list(Winf,K,t0),givenMeanWeight=NULL){

                              #- Make object to store new weights in
                              storeWt           <- array(NA,dim=c(length(ages),length(projPeriod),1,1,1,iters))
                              #- Get historic observed weights
                              wtsdat            <- as.data.frame(histWeight[,histPeriod])
                              colnames(wtsdat)  <- c(colnames(wtsdat)[1:6],"weight")
                              #- Fit von Bertalanffy to data
                              wvBertalanffy     <- nls(weight ~   (Winf*(1-exp(-K*(age-t0)))^3.038),data=wtsdat,start=list(Winf=estim$Winf,K=estim$K,t0=estim$t0),control=list(minFactor=1/2^16,maxiter=1000,warnOnly=T))

                              #- Get mean and sd of fit
                              meanAge <- log(predict(wvBertalanffy,new=data.frame(age=ages)))
                              if(is.null(givenMeanWeight)==F) meanAge <- log(givenMeanWeight)
                              sdAge   <- apply(log(histWeight),1,sd,na.rm=T)@.Data
                              #- Determine 2sd bound
                              bounds  <- cbind(exp(meanAge)-2*apply(histWeight,1,sd)@.Data,
                                               exp(meanAge)+2*apply(histWeight,1,sd)@.Data)
                              bounds[bounds<=0] <- 0.01

                              #- Draw new weights
                              for(iAge in 1:length(ages)){
                                for(iYr in 1:length(projPeriod)){
                                  for(iTer in 1:iters){
                                    set.seed(iAge*iYr*iTer)
                                    if(iAge == 1){ storeWt[iAge,iYr,,,,iTer] <- rtlnorm(1,mean=meanAge[iAge],sd=sdAge[iAge,,,,,],lower=bounds[iAge,1],upper=bounds[iAge,2])
                                    } else {       storeWt[iAge,iYr,,,,iTer] <- rtlnorm(1,mean=meanAge[iAge],sd=sdAge[iAge,,,,,],lower=storeWt[(iAge-1),iYr,,,,iTer],upper=bounds[iAge,2])
                                      }
                                  }
                                }
                              }
                            return(storeWt)}
                            
#-------------------------------------------------------------------------------
#- Update stock objects with 1 year of fisheries and biological data
#-------------------------------------------------------------------------------
updateStocks <- function(iStocks,iFishery,yr,iBiol,iCatchResids){


# iStocks        <- tmp_stocks
# iFishery       <- tmp_fishery
# yr             <- yrmin1
# iBiol          <- tmp_biol
# iCatchResids   <- ctch 
#



                  iStocks@catch.n[,ac(yr)]        <- unitSums(catch.n(iFishery)[,ac(yr)])   * iCatchResids[,ac(yr)]
                  iStocks@landings.n[,ac(yr)]     <- unitSums(landings.n(iFishery)[,ac(yr)])* iCatchResids[,ac(yr)]
                  iStocks@discards.n[,ac(yr)]     <- unitSums(discards.n(iFishery)[,ac(yr)])
                  iStocks@landings.wt[,ac(yr)]    <- unitSums(landings.n(iFishery)[,ac(yr)] * landings.wt(iFishery)[,ac(yr)]) / unitSums(landings.n(iFishery)[,ac(yr)])
                  iStocks@discards.wt[,ac(yr)]    <- 0
                  iStocks@catch.wt[,ac(yr)]       <- iStocks@landings.wt[,ac(yr)]
                  iStocks@catch[,ac(yr)]          <- computeCatch(iStocks)[,ac(yr)]
                  iStocks@landings[,ac(yr)]       <- computeLandings(iStocks)[,ac(yr)]
                  iStocks@discards[,ac(yr)]       <- computeDiscards(iStocks)[,ac(yr)]
                  
                  iStocks@mat[,ac(yr)]            <- iBiol@fec[,ac(yr)]
                  iStocks@stock.wt[,ac(yr)]       <- iBiol@wt[,ac(yr)]
                  iStocks@m[,ac(yr)]              <- iBiol@m[,ac(yr)]
                return(iStocks)}
                
#-------------------------------------------------------------------------------
#- Calculate the harvest by fleet given a TAC for a stock
#-------------------------------------------------------------------------------
fleet.harvest <- function(stk,iYr,TACS){
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter) res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleF,stk=iter(stk,iTer),iYr=iYr,TACS=c(iter(TACS,iTer)),nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL)$par
                    stk@harvest[,iYr]               <- sweep(stk@harvest[,iYr],c(3,6),res,"*")
                 return(stk@harvest[,iYr])}

  #faster version
fleet.harvestF<- function(stk,iYr,TACS){
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter){
                      Ns      = stk@stock.n[,iYr,1,,,iTer]@.Data
                      Fs      = stk@harvest[,iYr,,,,iTer]@.Data
                      Cwts    = stk@catch.wt[,iYr,,,,iTer]@.Data
                      Ms      = stk@m[,iYr,1,,,iTer]@.Data
                      res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=TACS[,,,,,iTer],nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL)$par
                    }
                    stk@harvest[,iYr]               <- sweep(stk@harvest[,iYr],c(3,6),res,"*")
                 return(stk@harvest[,iYr])}
                 
  #faster version
fleet.harvestFF<- function(stk,iYr,TACS){
                    lin     <- substr(R.Version()$os,1,3)== "lin"
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter){
                      Ns      = c(stk@stock.n[,iYr,1,,,iTer]@.Data)
                      Fs      = stk@harvest[,iYr,,,,iTer,drop=T]@.Data
                      Cwts    = stk@catch.wt[,iYr,,,,iTer,drop=T]@.Data
                      Ms      = c(stk@m[,iYr,1,,,iTer]@.Data)
                      if(lin) res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=c(TACS[,,,,,iTer]),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits),nls.lm.control(ftol = (.Machine$double.eps)))$par
                      if(!lin)res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=c(TACS[,,,,,iTer]),jac=NULL,nls.lm.control(ftol = (.Machine$double.eps)))$par
                    }
                    stk@harvest[,iYr]               <- sweep(stk@harvest[,iYr],c(3,6),res,"*")
                 return(stk@harvest[,iYr])}
                 
#-------------------------------------------------------------------------------
#- Function to scale the F pattern for stock
#-------------------------------------------------------------------------------
rescaleF      <- function(mult,stk.=stk,iYr.=iYr,TACS.=TACS){
                    stk.@harvest[,iYr.] <- sweep(stk.@harvest[,iYr.],3,mult,"*")
                    stkZ                <- unitSums(stk.@harvest[,iYr.]) + stk.@m[,iYr.,1]
                    res                 <- sqrt(c((TACS. - c(apply(sweep(stk.@stock.n[,iYr.] * stk.@catch.wt[,iYr.] * sweep(stk.@harvest[,iYr.],c(1:2,4:6),stkZ,"/"),c(1:2,4:6),(1-exp(-stkZ)),"*"),3:6,sum,na.rm=T)))^2))
                 return(res)}

  #Faster version
rescaleFF     <- function(mult,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS.=TACS){
                    Fs    <- sweep(Fs,3,mult,"*")
                    stkZ  <- apply(Fs,1,sum) + Ms
                    res   <- sqrt(c((TACS. - c(apply(sweep(sweep(Fs,c(1:2,4:6),stkZ,"/")* Cwts,c(1:2,4:6),Ns * (1-exp(-stkZ)),"*") ,3:6,sum))))^2)
                 return(res)}
                 
  #Faster version
rescaleFFF     <- function(mult,Ns.=Ns,Fs.=Fs,Cwts.=Cwts,Ms.=Ms,TACS.=TACS){
                    Fs.    <- t(t(Fs.) * mult)
                    stkZ  <- rowSums(Fs.) + Ms.
                    res   <- sqrt(c((TACS. - c(colSums(sweep(sweep(Fs.,1,stkZ,"/")* Cwts.,1,Ns. * (1-exp(-stkZ)),"*")))))^2)
                 return(res)}
                 
#-------------------------------------------------------------------------------
#- Calculate the catch by fleet
#-------------------------------------------------------------------------------
harvestCatch  <-  function(stk.,iYr){
                    stkZ      <- unitSums(stk.@harvest[,iYr]) + stk.@m[,iYr,1]
                    res       <- apply(sweep(stk.@stock.n[,iYr] * stk.@catch.wt[,iYr] * sweep(stk.@harvest[,iYr],c(1:2,4:6),stkZ,"/"),c(1:2,4:6),(1-exp(-stkZ)),"*"),3:6,sum,na.rm=T)
                  return(res)}
                  
#-------------------------------------------------------------------------------
#- Management plan: calculate A and B TAC
#-------------------------------------------------------------------------------

                 
find.FABF     <- function(mult,stk.,mpPoints.){

#mult<-1
#stk.<- iter(stf[,FcY],iTer)
#mpPoints.<-mpPoints
                    Fs      <- sweep(stk.@harvest@.Data,3,mult,"*")
                    Ns      <- stk.@stock.n@.Data[,,1,,,]
                    Ms      <- stk.@m@.Data[,,1,,,]
                    Hspwns  <- stk.@harvest.spwn@.Data[,,1,,,]
                    Mspwns  <- stk.@m.spwn@.Data[,,1,,,]
                    Mats    <- stk.@mat@.Data[,,1,,,]
                    Swghts  <- stk.@stock.wt@.Data[,,1,,,]
                    Cwghts  <- stk.@catch.wt@.Data

                    bigF              <- apply(Fs,1,sum)
                    ssb               <- sum(Ns * Swghts * exp(-bigF*Hspwns - Ms*Mspwns) * Mats)
                    if(ssb < mpPoints.$Btrigger) resA <- mpPoints.$Ftarget* ssb / mpPoints.$Btrigger
                    if(ssb > mpPoints.$Btrigger) resA <- mpPoints.$Ftarget   
                                 
                    fbarA     <- mean(bigF[ac(4:8)])
                    
                    ret       <- (fbarA-resA)^2
               
                 return(ret)}

#-------------------------------------------------------------------------------
#- Function to calculate the change in effort
#-------------------------------------------------------------------------------
f21t <- function(fmult,Q1,morts,Fsq,landsel,nums,land_wts){
    bigF <- Fsq * fmult
    return(abs(Q1 - sum( ((bigF)/(bigF + morts))*nums*(1-exp(-bigF-morts))* land_wts)))
  }

#-------------------------------------------------------------------------------
#- Function to return the change in effort
#-------------------------------------------------------------------------------
f31t <- function(Q1,iBiol,ImY,iFishery,iFailRuns=NULL){
  lan_wt  <- landings.wt( iFishery)[,ac(ImY)]
  fpattern<- catch.sel(iFishery)[,ac(ImY)]
  iters   <- dims(Q1)$iter
  nUnits  <- dims(Q1)$unit
  fm      <- matrix(NA,nrow=dims(iFishery)$unit,ncol=iters,dimnames=list(dimnames(iFishery@landings.n)$unit,dimnames(iFishery@landings.n)$iter))
  if(is.null(iFailRuns)==T){
    for(i in 1:iters){
      iMort    <- iBiol@m[,ac(ImY),,,,i]@.Data
      iNum     <- iBiol@n[,ac(ImY),,,,i]@.Data
      iFsq     <- fpattern[,,,,,i]@.Data
      iLanwt   <- lan_wt[,,,,,i]@.Data
      iQ1      <- Q1[,ac(ImY),,,,i]@.Data
      fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiol,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
    }
  } else {
    for(i in (1:iters)[-iFailRuns]){
      iMort    <- iBiol@m[,ac(ImY),,,,i]@.Data
      iNum     <- iBiol@n[,ac(ImY),,,,i]@.Data
      iFsq     <- fpattern[,,,,,i]@.Data
      iLanwt   <- lan_wt[,,,,,i]@.Data
      iQ1      <- Q1[,ac(ImY),,,,i]@.Data
      fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiol,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
    }
  }
  return(fm)
}






f31tF <- function(Q1,iBiol,ImY,iFishery,iFailRuns=NULL){

 #  Q1=TAC*TACusage ; iBiol=biol ; iFishery =fishery


  lin     <- substr(R.Version()$os,1,3)== "lin"
  lan_wt  <- landings.wt( iFishery)[,ac(ImY)]
  fpattern<- catch.sel(iFishery)[,ac(ImY)]
  iters   <- dims(Q1)$iter
  nUnits  <- dims(Q1)$unit
  fm      <- matrix(NA,nrow=dims(iFishery)$unit,ncol=iters,dimnames=list(dimnames(iFishery@landings.n)$unit,dimnames(iFishery@landings.n)$iter))
  if(is.null(iFailRuns)==T){
    for(i in 1:iters){
      iMort    <- c(iBiol@m[,ac(ImY),,,,i]@.Data)
      iNum     <- c(iBiol@n[,ac(ImY),,,,i]@.Data)
      iFsq     <- fpattern[,,,,,i,drop=T]@.Data
      iLanwt   <- lan_wt[,,,,,i,drop=T]@.Data
      iQ1      <- c(Q1[,ac(ImY),,,,i]@.Data)
      if(lin) fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
      if(!lin)fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL)$par
    }
  } else {
    for(i in (1:iters)[-iFailRuns]){
      iMort    <- c(iBiol@m[,ac(ImY),,,,i]@.Data)
      iNum     <- c(iBiol@n[,ac(ImY),,,,i]@.Data)
      iFsq     <- fpattern[,,,,,i,drop=T]@.Data
      iLanwt   <- lan_wt[,,,,,i,drop=T]@.Data
      iQ1      <- c(Q1[,ac(ImY),,,,i]@.Data)
      if(lin) fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL,lower=rep(0,nUnits),upper=rep(1e5,nUnits))$par
      if(!lin)fm[,i]   <- nls.lm(par=rep(1,nUnits),rescaleFBiolF,iFsq.=iFsq,iMort.=iMort,iNum.=iNum,iLanwt.=iLanwt,TACS.=iQ1,nls.lm.control(ftol = (.Machine$double.eps)),jac=NULL)$par
    }
  }
  return(fm)
}


#-------------------------------------------------------------------------------
#- Function to scale the F pattern for biol
#-------------------------------------------------------------------------------
rescaleFBiol <- function(mult,iFsq.,iMort.,iNum.,iLanwt.,TACS.){
                    iF          <- sweep(iFsq.,3,mult,"*")
                    stkZ        <- apply(iF,1,sum) + iMort.
                    res         <- sqrt(c((c(TACS.) - c(apply(sweep(sweep(iLanwt. * iF,c(1:2,4:6),stkZ,"/"),c(1:2,4:6),iNum. * (1-exp(-stkZ)),"*"),3:6,sum,na.rm=T)))^2))
                 return(res)}

rescaleFBiolF <- function(mult,iFsq.,iMort.,iNum.,iLanwt.,TACS.){
                    iF          <- t(t(iFsq.) * mult)
                    stkZ        <- rowSums(iF) + iMort.
                    res         <- sqrt(c((TACS. -    c(colSums(sweep(sweep(iF,1,stkZ,"/")* iLanwt.,1,iNum. * (1-exp(-stkZ)),"*")))))^2)
                 return(res)}


#-------------------------------------------------------------------------------
#- Function to calculate the ssb based on biological numbers at age
#-------------------------------------------------------------------------------

ssbb <- function(iBiol,iHarvest,iStck){

#  iBiol= biol[,ac(yrs)] ; iHarvest =   ;iStck=
          res <- quantSums(iBiol@n * iBiol@wt * exp(-iHarvest * iStck@harvest.spwn - iBiol@m * iStck@m.spwn) * iBiol@fec)
        return(res)}

#-------------------------------------------------------------------------------
#- Function to estimate the mean and sd of the lognormal distribution that fits
#  the observed recruitment best
#-------------------------------------------------------------------------------

optimRecDistri <- function(pars,recs){
                    sim <- plnorm(recs,meanlog=pars[1],sdlog=pars[2])
                    lik <- sum(abs(log(sim) - log(seq(0,1,length.out=length(recs)+2)[-c(1,length(recs)+2)])))
                  return(lik)}

#-------------------------------------------------------------------------------
#- Function for paralelisation
#-------------------------------------------------------------------------------

startCluster <- function(snow.cores){
                  #-Start up cluster
                  library(snow)
                  library(Rmpi)
                  num.cores       <- snow.cores
                  snow.cluster    <- c(rep("localhost",num.cores))
                  cl              <- makeMPIcluster(count=length(snow.cluster))
                  options(timeout=14*86400)
                  clusterEvalQ(cl,library(FLSAM))
                  clusterEvalQ(cl,library(MASS))
                  clusterEvalQ(cl,library(msm))
                return(cl)}
                
exportCluster <- function(cl,lst){
                  clusterExport(cl,lst)
                 }

retroLB                     <- function(x){
                                cat(paste("\n\n\n\ THIS IS ITERATION NUMBER",x,"\n\n\n"))
                                stck  <- window(iter(stocks,x),end=histMaxYr)
                                tun   <- FLIndices()
                                for(iSurv in names(surveys))
                                  tun[[iSurv]] <- iter(window(surveys[[iSurv]],end=range(NSH.tun[[iSurv]])["maxyear"]),x)
                                res <- retro(stock=stck,indices=tun,
                                             control=ctrl,retro=10,base.assess=base.assess)
                               return(res)}

#-------------------------------------------------------------------------------
#- Function to overwrite the updated tmp_stocks to stocks (but do not overwrite
#   the elements listed below from tmp to stocks, hence they are taken from
#   the stocks object)
#-------------------------------------------------------------------------------


tmp2stocks <- function(stocks,tmp_stocks,yr){
                TaYtmp_stocks         <- tmp_stocks[,ac(yr)]
                TaYstocks             <- stocks[,ac(yr)]
                TaYtmp_stocks@stock   <- TaYstocks@stock
                TaYtmp_stocks@stock.n <- TaYstocks@stock.n
                TaYtmp_stocks@harvest <- TaYstocks@harvest
                TaYtmp_stocks@name    <- TaYstocks@name
                TaYtmp_stocks@desc    <- TaYstocks@desc
                TaYtmp_stocks@range   <- TaYstocks@range
                stocks[,ac(yr)]       <- TaYtmp_stocks
              return(stocks)}
              
              
              
#-----------------------------------------------------------------------------------------------------------
# Recruitment function => gives the R accoring to the SSB and SR parameters
#-----------------------------------------------------------------------------------------------------------

B2R<-function (Biom, SRp,ite,yr,recreg){

#Biom<- ssb(biol)[,ac(iYr-1),,,,its]
#SRp<-  SRmod[its,]
#ite<-  its
#yr<-   iYr
#recreg<- RecType
#
      if( recreg  !=  "SGRec"){
          mod <-  SRp$mod
          A     <-  SRp$A
          B     <-  SRp$B
          sig   <-  SRp$sigma  # this rescaling should have been done in the Bayesian protocol

          if (mod=="segreg" ) mu <-  A*(Biom>=B)+A*Biom*(Biom<B)/B
          if (mod=="ricker" ) mu <-  A*Biom*exp(-B*Biom)
          if (mod=="bevholt") mu <-  A*Biom/(B+Biom)
          
          REC<-exp(log(mu)+rnorm(1,0,sig)) 
          
          REC<- max(REC,10)
        }
        if ( recreg == "SGRec")
          REC   <- SRp[1,yr]

    return(REC)}              