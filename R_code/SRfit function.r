###########################################################################################################################
# fitting with lognormal errors
###########################################################################################################################



dats<-data.frame(R=R,SSB=SSB)

### seg reg model
#------------------------------
#try a grid of starting values
LL<-c()
a <-mean(dats$R)
b <- mean(dats$SSB)
sigma <-sd(log(dats$R))
ini<-expand.grid(seq(as.numeric(a)*0.5,as.numeric(a)*1.5,length.out=5),seq(as.numeric(b)*0.5,as.numeric(b)*1.5,length.out=5),sigma*10^seq(-1,1,length.out=3))
rm(a,b,sigma)
LL<-c()
mods<-list()
for( i in 1:dim(ini)[1])
{
Ll<-NA
mod<-"pas de solution"
try(mod<-optim(ini[i,],fn=logL.logN,SRdats=dats,mod="segreg",gr=NULL,lower=c(0,0,0), method= "Nelder-Mead"),silent=T)
try(Ll<-mod$value,silent=T)
mods[[i]]<-mod
LL<-c(LL,Ll)
}
ini<-cbind(ini,LL)
ini<-ini[!is.na(ini$LL),]
best<-ini[ini$LL== min(ini$LL),]
best<-best[1,]

#do the final run using the best starting values
fitSG.logN<-optim(best[1:3],fn=logL.logN,SRdats=dats,mod="segreg",lower=c(0,0,0),gr=NULL, method= "Nelder-Mead")

rm(ini,best,LL)



### Ricker model
#------------------------------
#try a grid of starting values

ini<-lm(log(R/SSB)~SSB,data=dats)
a <-exp(ini$coefficients[1])
b <- - ini$coefficients[2]
sigma <- sd(log(dats$R))
ini<-expand.grid(seq(as.numeric(a)*0.5,as.numeric(a)*1.5,length.out=5),seq(as.numeric(b)*0.5,as.numeric(b)*1.5,length.out=5),sigma*10^seq(-2,2,length.out=3))
rm(a,b,sigma)
LL<-c()
mods<-list()
for( i in 1:dim(ini)[1])
{
Ll<-NA
mod<-"pas de solution"
try(mod<-optim(ini[i,],fn=logL.logN,SRdats=dats,mod="ricker",gr=NULL, method= "Nelder-Mead"),silent=T)
try(Ll<-mod$value,silent=T)
mods[[i]]<-mod
LL<-c(LL,Ll)
}
ini<-cbind(ini,LL)
ini<-ini[!is.na(ini$LL),]
best<-ini[ini$LL== min(ini$LL),]
best<-best[1,]


#do the final run using the best starting values
fitRi.logN<-optim(best[1:3],fn=logL.logN,SRdats=dats,mod="ricker",gr=NULL, method= "Nelder-Mead")


### bevholt model
#------------------------------
#try a grid of starting values
LL<-c()
invR<-1/dats$R
invS<-1/dats$SSB
ini<-lm(invR~invS)
a <-1/(ini$coefficients[1])
b <-  ini$coefficients[2] *a
sigma <-sd(log(dats$R))
ini<-expand.grid(seq(as.numeric(a)*0.5,as.numeric(a)*1.5,length.out=5),seq(as.numeric(b)*0.5,as.numeric(b)*1.5,length.out=5),sigma*10^seq(-2,2,length.out=3))
rm(a,b,sigma)
LL<-c()
mods<-list()
for( i in 1:dim(ini)[1])
{
Ll<-NA
mod<-"pas de solution"
try(mod<-optim(ini[i,],fn=logL.logN,SRdats=dats,mod="bevholt",gr=NULL,lower=c(0,0,0), method= "Nelder-Mead"),silent=T)
try(Ll<-mod$value,silent=T)
mods[[i]]<-mod
LL<-c(LL,Ll)
}
ini<-cbind(ini,LL)
ini<-ini[!is.na(ini$LL),]
best<-ini[ini$LL== min(ini$LL),]
best<-best[1,]


#do the final run using the best starting values
try(fitBH.logN<-optim(best[1:3],fn=logL.logN,SRdats=dats,mod="bevholt",lower=c(0,0,0),gr=NULL, method= "Nelder-Mead"),silent=T)

rm(ini,best,LL)




# combine results and find out which SR form is best
rm(res)

res<-data.frame(model=NA,distri=NA,a=NA,b=NA,sigma=NA,lLik=NA,AICc=NA)
try(res<-rbind(res,c("Ricker","logNormal",fitRi.logN$par[1],fitRi.logN$par[2],fitRi.logN$par[3],-fitRi.logN$value,NA)),silent=T)
try(res<-rbind(res,c("BevHolt","logNormal",fitBH.logN$par[1],fitBH.logN$par[2],fitBH.logN$par[3],-fitBH.logN$value,NA)),silent=T)
try(res<-rbind(res,c("Hockey Stick","logNormal",fitSG.logN$par[1],fitSG.logN$par[2],fitSG.logN$par[3],-fitSG.logN$value,NA)),silent=T)
res<-res[-1,]
res$a<-as.numeric(res$a);res$b<-as.numeric(res$b);res$sigma<-as.numeric(res$sigma);res$lLik<-as.numeric(res$lLik);
k<-3
res$AICc<-2*k-2*res$lLik+2*k *(k+1)/(dim(dats)[1]-k-1)

res<-res[res$AIC==min(res$AIC),]
rm(fitRi.logN)
rm(fitBH.logN)
rm(fitSG.logN)

