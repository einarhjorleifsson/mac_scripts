idxssb <- which(HCRSSB[,ac(2013)] < mpPoints$Btrigger)
(mpPoints$FadultA - mpPoints$FadultB) / ((mpPoints$Btrigger -mpPoints$Blim) / 1e6)*((HCRSSB[,ac(2013),,,,idxssb]-mpPoints$Blim)/1e6) + mpPoints$FadultB


iTAC <- TAC
idxssb <- which(HCRSSB[,ac(2013)] > mpPoints$Btrigger)
idxtac <- which(HCRTAC[,ac(2013),c("A")] > 1.15*iTAC[,ac(2012),c("A")])
length(unique(c(idxssb,idxtac)))

length(which(HCRTAC[,ac(2013),c("A")] > 1.15*iTAC[,ac(2012),c("A")]))


by Btrigger = 1.7 is de TAC voorgesteld een beetje lager dan bij Btrigger = 1.6.
Maar doordat er dezelfde startwaarde is, kan de F van de B vloot iets hoger zijn
bij een Btrigger van 1.7.

Van de 831 die boven Btrigger=1.7 zitten, zijn er 453 die ook tegen de 15% IAV TAC aanlopen.
Van de 910 die boven Btrigger=1.6 zitten, zijn er 454 die ook tegen de 15% IAV TAC aanlopen.

Hoeveel SSB > Btrigger en TAC > 1.15*TAC

length(which(HCRSSB[,ac(2013)] > mpPoints$Btrigger & HCRTAC[,ac(2013),c("A")] > 1.15*iTAC[,ac(2012),c("A")]))
Dat zijn er altijd minder dan 500. Oftewel, wat je aan unconstraint TAC krijgt is belangrijker.
En ongeacht je lager gaat in Btrigger neemt het aantal runs dat constraint is by IAV niet toe.


which(round(apply(unitSums(stf@harvest[f26,ac(2013)]),2:6,mean),5) < 0.24 &
      )
      
      
#1900000
idxssb <- which(HCRSSB[,ac(2013)] > mpPoints$Btrigger) #571 zijn groter dan Btrigger
c(apply(unitSums(stf@harvest[f26,ac(2013),,,,idxssb]),2:6,mean)) #ze zijn allemaal 0.24

#van idxssb zijn er 423 die tegen TAC limit aanlopen
length(which(HCRTAC[,ac(2013),1,,,idxssb] > 1.15*405000))
#van idxssb zijn er 148 die niet tegen TAC limit aanlopen
#oftewel: 1000-423 = 577 worden door F bepaald, niet door limit

#1300000
idxssb <- which(HCRSSB[,ac(2013)] > mpPoints$Btrigger) #998 zijn groter dan Btrigger
c(apply(unitSums(stf@harvest[f26,ac(2013),,,,idxssb]),2:6,mean)) #ze zijn niet allemaal 0.24
#van idxssb zijn er 454 die tegen TAC limit aanlopen
length(which(HCRTAC[,ac(2013),1,,,idxssb] > 1.15*405000))


plot(x=HCRSSB[,ac(2013)],y=quantMeans(unitSums(stf@harvest[f26,ac(2013)])))




#
idxbtrigger <- which(HCRSSB > mpPoints$Btrigger)
which(HCRTAC[,ac(2013),1,,,idxbtrigger] < 1.15*405000)


#1.9e6
idx <- which(HCRTAC[,ac(2013),1] %in% sort(c(HCRTAC[,ac(2013),1]))[500:501])
print(idx)
c(HCRSSB[,ac(2013),,,,idx])
c(apply(unitSums(stf@harvest[f26,FcY,,,,idx]),2:6,mean))
c(HCRTAC[,ac(2013),1,,,idx])


