################################################################################
# NSH_SAM Control for Assessment
#
# $Rev: 705 $
# $Date: 2012-02-14 19:02:57 +0100 (di, 14 feb 2012) $
#
# Author: HAWG model devlopment group
#
# Sets up a control object for use by Step 04 assessments i.e. the "refined data" run
#
# Developed with:
#   - R version 2.13.0
#   - FLCore 2.4
#
# To be done:
#
# Notes: Have fun running this assessment!
#
################################################################################

Mac.ctrl <- FLSAM.control(Mac,Mac.tun)

#Set the variances. Separate variance for recruitment and plus group
#Fishing mortality RWs are set from an analysis of ICA VPA results
Mac.ctrl@logN.vars[]      <- c(1,rep(2,dims(Mac)$age-1))
Mac.ctrl@f.vars["catch",] <- 1

#All fishing mortality states are free except
#oldest ages to ensure stablity
Mac.ctrl@states["catch",] <- seq(dims(Mac)$age)
Mac.ctrl@states["catch",ac(7:12)] <- 101

#Group observation variances of catches to ensure stability
Mac.ctrl@obs.vars["catch",]  <- 1

Mac.ctrl@obs.vars["R-idx(log transf)",ac(0)]    <- 2
Mac.ctrl@obs.vars["Swept-idx",ac(6:11)]  <- 3

Mac.ctrl@catchabilities["Swept-idx",ac(6:11)]<-2

#Finalise
Mac.ctrl@name <- "Final Assessment"
Mac.ctrl <- update(Mac.ctrl)

