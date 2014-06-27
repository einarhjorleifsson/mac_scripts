codePath      <- "N:/Projecten/ICES WG/WKHELP/R/"
if(substr(R.Version()$os,1,3)== "lin")
  codePath    <- sub("N:/","/media/n/",codePath)

Fmsy     <- 0.3
mpPoints <- list(Blim           =0.8e6,
                 Bpa            =1.3e6,
                 Btrigger       =1.5e6,
                 FadultA        =Fmsy,
                 FadultB        =0.1,
                 FjuvA          =0.05,
                 FjuvB          =0.04,
                 stabilityBreak ="Blim",
                 scen           ="LTMP")
                 
#-------------------------------------------------------------------------------
#- Banking
#-------------------------------------------------------------------------------
mpPoints$scen <- "Bank"
scen <- mpPoints$scen
for(iBreak in seq(1.3e6,2.0e6,0.1e6)){
  for(iFmsy in c(0.24,0.3)){
    mpPoints$Btrigger <- iBreak
    mpPoints$FadultA  <- iFmsy
    source(file.path(codePath,"03b_runMSE.r"))
  }
}#10 * 1uur
  
#-------------------------------------------------------------------------------
#- BB
#-------------------------------------------------------------------------------

mpPoints$scen <- "BB"
scen <- mpPoints$scen
for(iBreak in seq(1.3e6,2.0e6,0.1e6)){
  for(iFmsy in c(0.24,0.3)){
    mpPoints$Btrigger <- iBreak
    mpPoints$FadultA  <- iFmsy
    source(file.path(codePath,"03b_runMSE.r"))
  }
}# 10 * 1uur

#-------------------------------------------------------------------------------
#- F3years
#-------------------------------------------------------------------------------
mpPoints$scen <- "F3years"
scen <- mpPoints$scen
for(iStable in c("Blim","Btrigger")){
  for(iFmsy in c(0.24,0.3)){
    mpPoints$FadultA        <- iFmsy
    mpPoints$stabilityBreak <- iStable
    source(file.path(codePath,"03b_runMSE.r"))
  }
}# 4 * 1uur

#-------------------------------------------------------------------------------
#- Bycatch
#-------------------------------------------------------------------------------
mpPoints$scen   <- "LTMP"
scen <- "juv"
for(iJuv in c(1,0.75)){
  for(iAdu in c(0.2,0.24,0.3)){
    if(iAdu == 0.2 & iJuv == 1){
      mpPoints$FadultA  <- iAdu
      mpPoints$FjuvA    <- 0.05 *iJuv
      mpPoints$FjuvB    <- 0.04 *iJuv
    }
    if(iAdu %in% c(0.24,0.3) & iJuv == 0.75){
      mpPoints$FadultA  <- iAdu
      mpPoints$FjuvA    <- 0.05 *iJuv
      mpPoints$FjuvB    <- 0.04 *iJuv
    }
    source(file.path(codePath,"03b_runMSE.r"))
  }
}# 3 * 1uur

#-------------------------------------------------------------------------------
#- LTMP
#-------------------------------------------------------------------------------
mpPoints$scen   <- "LTMP"
scen <- mpPoints$scen
for(iAdu in c(0.2,0.24,0.25,0.3)){
  mpPoints$FadultA  <- iAdu
  for(iStable in c("Blim","Btrigger")){
    mpPoints$stabilityBreak <- iStable
    if(iStable == "Blim"){
      if(iAdu %in% c(0.24,0.3)){
        for(iBreak in seq(1.3e6,2.0e6,0.1e6)){
          mpPoints$Btrigger <- iBreak
          source(file.path(codePath,"03b_runMSE.r"))
        }
      } else {
          mpPoints$Btrigger <- 1.5e6
          source(file.path(codePath,"03b_runMSE.r"))
        }
    } else {
        mpPoints$Btrigger <- 1.5e6
        source(file.path(codePath,"03b_runMSE.r"))
      }
  }
}#Btrigger = 4, Blim=4 = 2+2*5: total 16 * 1uur

#-------------------------------------------------------------------------------
#- meanTAC
#-------------------------------------------------------------------------------
mpPoints$scen   <- "meanTAC"
scen <- mpPoints$scen
for(iAdu in c(0.24,0.3)){
  mpPoints$FadultA <- iAdu
  for(iStable in c("Blim","Btrigger")){
    mpPoints$stabilityBreak <- iStable
    if(iStable == "Blim"){
      for(iBreak in seq(1.3e6,2.0e6,0.1e6)){
        mpPoints$Btrigger <- iBreak
        source(file.path(codePath,"03b_runMSE.r"))
      }
    } else {
        mpPoints$Btrigger <- 1.5e6
        source(file.path(codePath,"03b_runMSE.r"))
      }
  }
}# 2*5 + 2: total 12 * 1uur

#-------------------------------------------------------------------------------
#- noIAV
#-------------------------------------------------------------------------------
mpPoints$scen   <- "noIAV"
scen <- mpPoints$scen
for(iAdu in c(0.24,0.3)){
  mpPoints$FadultA <- iAdu
  source(file.path(codePath,"03b_runMSE.r"))
}# 2 * 1uur

#-------------------------------------------------------------------------------
#- FIAV
#-------------------------------------------------------------------------------
mpPoints$scen   <- "FIAV"
scen <- mpPoints$scen
for(iAdu in c(0.24,0.3)){
  mpPoints$FadultA <- iAdu
  for(iStable in c("Blim","Btrigger")){
    mpPoints$stabilityBreak <- iStable
    if(iStable == "Blim"){
      for(iBreak in seq(1.3e6,2.0e6,0.1e6)){
        mpPoints$Btrigger <- iBreak
        source(file.path(codePath,"03b_runMSE.r"))
      }
    } else {
        mpPoints$Btrigger <- 1.5e6
        source(file.path(codePath,"03b_runMSE.r"))
      }
  }
}# 2*5+2: total 12 * 1uur



