#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 05-Jun-2012
#
# Build for R2.13.2, 32bits
#-------------------------------------------------------------------------------

#- Scenario descriptions
Fmsy          <- 0.25
                      
LTMP          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.36e6,
                                Ftarget=0.2,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.36e6,
                                Ftarget=0.3,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.36e6,
                                Ftarget=0.4,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.36e6,
                                Ftarget=0.25,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.36e6,
                                Ftarget=0.35,
                                stabilityBreak="Btrigger",scen="LTMP"))
                                
                                
                                
              