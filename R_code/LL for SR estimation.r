# define the loglikelihood functions


logL.N<-function(pars,SRdats,mod)
{
ai<-pars[1] ;bi<-pars[2];sigmai<-pars[3]
ssb<-SRdats$SSB
r<-SRdats$R
if(mod=="ricker")  rmod <- ai*ssb*exp(-bi*ssb)
if(mod=="bevholt") rmod <- ai*ssb /(bi + ssb)
if(mod=="segreg")
    {
    rmod<-rep(NA,length(r))
    rmod[ssb<bi]<-ai * ssb [ssb<bi] / bi
    rmod[ssb>=bi]<- ai
    }
like<-dnorm(r,rmod,rmod*sigmai)      # I have checked and this is the same as if writen explicitely the pdf of the normal distribtuion
loglike<-sum(log(like))
return(- loglike)
}



logL.N2<-function(pars,SRdats,mod)
{
ai<-pars[1] ;bi<-pars[2];sigmai<-pars[3]
ssb<-SRdats$SSB
r<-SRdats$R
n<-length(r)
if(mod=="ricker")  rmod <- ai*ssb*exp(-bi*ssb)
if(mod=="bevholt") rmod <- ai*ssb /(bi + ssb)
if(mod=="segreg")
    {
    rmod<-rep(NA,length(r))
    rmod[ssb<bi]<-ai * ssb [ssb<bi] / bi
    rmod[ssb>=bi]<- ai
    }
like<-dnorm(r,rmod,rmod*sigmai)      # I have checked and this is the same as if writen explicitely the pdf of the normal distribtuion
loglike2<-sum(log(like))

like<-1/(sigmai*rmod*(2*pi)^0.5)*exp(-((r-rmod)^2)/(2*rmod^2*sigmai^2))
loglike<-sum(log(like))
cat(loglike,loglike2,"\n")
return(- loglike)

}




logL.logN<-function(pars,SRdats,mod)
{
ai<-pars[1];bi<-pars[2];sigmai<-pars[3]
ssb<-SRdats$SSB
r<-SRdats$R
n<-length(r)
if(mod=="ricker")  rmod <- ai*ssb*exp(-bi*ssb)
if(mod=="bevholt") rmod <- ai*ssb /(bi + ssb)
if(mod=="segreg")
    {
    rmod<-rep(NA,length(r))
    rmod[ssb<bi]<-  ai * ssb [ssb<bi] / bi
    rmod[ssb>=bi]<- ai
    }



like<-1/(r*sigmai*(2*pi)^0.5)*exp(-0.5*((log(r)-log(rmod))^2)/sigmai^2)
loglike<-sum(log(like))
return(- loglike)
}

logL.G<-function(pars,SRdats,mod)
{
ai<-pars[1] ;bi<-pars[2];scale<-pars[3]
ssb<-SRdats$SSB
r<-SRdats$R
if(mod=="ricker")  rmod <- ai*ssb*exp(-bi*ssb)
if(mod=="bevholt") rmod <- ai*ssb /(bi + ssb)
if(mod=="segreg")
    {
    rmod<-rep(NA,length(r))
    rmod[ssb<bi]<-ai * ssb [ssb<bi] / bi
    rmod[ssb>=bi]<- ai
    }

mu<-rmod
shape<-  mu/scale
like<-dgamma(r,shape=shape,scale=scale)
loglike<-sum(log(like))
return(- loglike)

}

