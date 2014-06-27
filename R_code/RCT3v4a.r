#usage: cod = RCT3("c:\\rct3\\data\\cod.csv",shrink=TRUE)
#usage: cod = RCT3("c:\\rct3\\data\\new.csv",old=FALSE,shrink=TRUE)

#filename = "c:\\rct3\\data\\PLE7AREC.csv"  #old style                                                                                                                    
#filename = "c:\\rct3\\data\\cod.csv"       #old style
#filename = "c:\\rct3\\data\\test.csv"      #old style
#filename = "c:\\rct3\\data\\new.txt"       #new style

#old=TRUE
#logged=NA
#minpoints=3
#regression = 'c'
#shrink = FALSE
#taper = 3
#tapert=20

RCT3 = function (filename, old=TRUE, logged = NA, minpoints = 3, regression='c', shrink=TRUE, taper=3, tapert=20, minse = 0, configfile = NA)
{
#RCT3v2.r
#
# TJE 19/01/2010
# Replicates the functionality of RCT3 under the following conditions:
#   - Negative survey or VPA data indicates no data
#   - No weightings applied to surveys
#   - Doesn't offer the option of excluding surveys greater than mean when not shrinking
#   - Needs testing with more general configfiles

# input parameters:
#  filename:     Path and filename of a file containing data as a csv including 'Nsurveys' surveys
#  old:          Boolean, TRUE (default) if file is old style, FALSE if new style
#  logged:       Boolean vector of length 'Nsurveys' (or a divisor of it) indicating whether survey data is already logged 
#                'NA' is equivalent to FALSE for all surveys
#  minpoints:    Minimum number of points required for a survey to be used
#  regression:   Character variable 'c' or 'p' indicating calibration or predictive regression respectively
#  shrink:       Whether shrinkage to mean is applied to the regression estimates. Generally (only) used with calibration regression
#  taper:        The order of the time weighted regression. Value must be one of 0, 1, 2 or 3
#  tapert:       Time that the taper occurs over
#  minse:        Minimum standard error to be applied to any survey  
#  configfile:   File where further parameters can be found, any other parameters given as arguments, except 'old'
#                are discarded, and values from configfile used instead

  #Reading in data 
  source(paste(codePath,"readrct3.r",sep="")) #For readRCT3, readConfigFile and sprintf function

  if (!is.na(configfile))
  {
    return(readConfigFile(filename, old, configfile))
  }

  .input = readRCT3(filename, old)
  if(is.null(.input)) return(NULL)
  IPhead = .input$head
  .input = .input$filecontents

  #Split up the data into calibration and prediction datasets
  .input[.input<0]=NA
  Calibration =  subset(.input, !is.na(VPA))
  Prediction = subset(.input, is.na(VPA))
  Nsurveys = dim(.input)[2]-2
  Ncalibration = dim(Calibration)[1]
  Nprediction = dim(Prediction)[1]    
  
  #Deal with error checking parameters
  if((is.na(sum(logged)))|(length(logged)==0))
  {
    logged = rep(FALSE, Nsurveys)
  } else {
    if (Nsurveys %% length(logged) != 0) warning("Nsurveys is not a multiple of length of argument 'logged'")
    logged = rep(logged, length.out=Nsurveys) 
  }
  
  for (s in 1:Nsurveys)
  {
    if (!logged[s])
    {
      Calibration[s+2] = log(Calibration[s+2]+1)
      Prediction[s+2] = log(Prediction[s+2]+1)
    }   
  }
  
    
  Ncalibration = dim(Calibration)[1]
  Nprediction = dim(Prediction)[1]   

  
 
  
  #####Check that minpoints in integer greater than or equal to 3
  #####Check that regression is 'p' or 'c', and this fits with shrink
  #####Check that taper is 0, 1, 2 or 3
  
  
  # Set up Output Matrices
  cn = c(colnames(Prediction)[1:Nsurveys+2],"mean")
  Gradient =  matrix(0,Nprediction,Nsurveys+1,dimnames=list(Prediction$Year,cn))
  Intercept = matrix(0,Nprediction,Nsurveys+1,dimnames=list(Prediction$Year,cn))
  StdErr =    matrix(0,Nprediction,Nsurveys,dimnames=list(Prediction$Year,colnames(Prediction)[1:Nsurveys+2]))
  Rsquare =   matrix(0,Nprediction,Nsurveys,dimnames=list(Prediction$Year,colnames(Prediction)[1:Nsurveys+2]))
  Npoints =   matrix(0,Nprediction,Nsurveys,dimnames=list(Prediction$Year,colnames(Prediction)[1:Nsurveys+2]))
  Index =     matrix(0,Nprediction,Nsurveys+1,dimnames=list(Prediction$Year,cn))
  Pred =      matrix(0,Nprediction,Nsurveys+1,dimnames=list(Prediction$Year,cn))
  StdErr2 =   matrix(0,Nprediction,Nsurveys+1,dimnames=list(Prediction$Year,cn))
  WAPW =      matrix(0,Nprediction,Nsurveys+1,dimnames=list(Prediction$Year,cn))
  .sumsqx =   matrix(0,Nprediction,Nsurveys,dimnames=list(Prediction$Year,colnames(Prediction)[1:Nsurveys+2]))
  xbar =      matrix(0,Nprediction,Nsurveys,dimnames=list(Prediction$Year,colnames(Prediction)[1:Nsurveys+2]))
  sumw =      matrix(0,Nprediction,Nsurveys,dimnames=list(Prediction$Year,colnames(Prediction)[1:Nsurveys+2]))  
  blank = rep(0,Nprediction)
  .YOutput = data.frame(Year=blank, WAPred=blank, LogPred=blank, ISE=blank, ESE=blank, VarRatio=blank)
  

#Fit linear model to survey and catch data
  for (year in 1:Nprediction)
  {
    for (s in 1:Nsurveys)
    {   
      if (!is.na(Prediction[year,s+2]))
      {
        x = Calibration[!is.na(Calibration[,s+2]),s+2]
        y = log(Calibration$VPA[!is.na(Calibration[,s+2])]+1)
        t = Calibration$Year[!is.na(Calibration[,s+2])]
        if (length(x) >= minpoints)
        {          
          w = (1-(((Prediction$Year[year]-t-1)/tapert)^taper))^taper    #Error if more than tapert years and taper = 0
          w[(Prediction$Year[year]-t-1)>tapert] = 0
          xbar[year,s] = weighted.mean(x,w)
          
          ybar = weighted.mean(y,w)
#          ybar = weighted.mean(log(Calibration$VPA[!is.na(.input[,2])]+1),(1-(((Prediction$Year[year]-Calibration$Year[!is.na(.input[,2])])/tapert)^taper))^taper)
          Sxx = sum(w*(x-xbar[year,s])^2)/sum(w)
          Sxy = sum(w*(x-xbar[year,s])*(y-ybar))/sum(w)
          Syy = sum(w*(y-ybar)^2)/sum(w)
          Gradient[year,s] = ifelse(regression=='c', Syy/Sxy, Sxy/Sxx)
          Intercept[year,s] = ybar - xbar[year,s]*Gradient[year,s]    
          StdErr[year,s] = sqrt(sum(w*(y - x*Gradient[year,s]-Intercept[year,s])^2)/(sum(w)-2)) 
          Rsquare[year,s] = Sxy*Sxy/(Sxx*Syy)
          Npoints[year,s] = length(x) 
    
          Index[year,s] = Prediction[year,s+2]
          Pred[year,s] = Index[year,s] * Gradient[year,s] + Intercept[year,s]
          .sumsqx[year,s] = sum(w*(x-xbar[year,s])^2) 
          sumw[year,s] = sum(w)
          StdErr2[year,s] = max(minse, StdErr[year,s]*sqrt(sum(w)/(sum(w)-2))*sqrt(1+(1/sum(w))+ (Index[year,s] - xbar[year,s])^2/.sumsqx[year,s]))     #STERR in the Fortran code
         
        } else {
          Gradient[year,s] = NA
          Intercept[year,s] = NA
          StdErr[year,s] = NA
          Rsquare[year,s] = NA
          Npoints[year,s] = NA
          Index[year,s] = NA
          Pred[year,s] = NA
          StdErr2[year,s] = NA 
        }
      } else {
          Gradient[year,s] = NA
          Intercept[year,s] = NA
          StdErr[year,s] = NA
          Rsquare[year,s] = NA
          Npoints[year,s] = NA
          Index[year,s] = NA
          Pred[year,s] = NA
          StdErr2[year,s] = NA 
      }


    } 
    w = (1-(((Prediction$Year[year]-Calibration$Year-1)/tapert)^taper))^taper 
    Gradient[year,Nsurveys+1] = 0 
    Intercept[year,Nsurveys+1] = weighted.mean(log(Calibration$VPA+1),w)
    Pred[year,Nsurveys+1] = weighted.mean(log(Calibration$VPA+1),w) 
    StdErr2[year,Nsurveys+1] = sqrt(sum(w*(log(Calibration$VPA+1)-Pred[year,Nsurveys+1])^2)/(sum(w)-1))
    Index[year,Nsurveys+1]=0
    
#stats for the overall prediction for the year    
    if (length(StdErr[year,!is.na(StdErr[year,])])>0)
    {
      if (shrink)
      {
        sumse2 = sum(StdErr2[year,]^-2,na.rm=TRUE) 
        WAPW[year,] = (StdErr2[year,]^-2)/sumse2
        .YOutput$LogPred[year] = sum(WAPW[year,]*Pred[year,],na.rm=TRUE)
        sumdiff = sum((Pred[year,]-.YOutput$LogPred[year])^2/(StdErr2[year,]^2),na.rm=TRUE)
        tmp = Prediction[year,-(1:2)]
        nsury = length(tmp[!is.na(tmp)]) + 1 #Number of surveys positively reporting in this year
        sinvsq = sum((StdErr2[year,]^-2),na.rm=TRUE) 
        .YOutput$ISE[year] = (sinvsq)^-0.5
        .YOutput$ESE[year] = sqrt(sumdiff / ((nsury-1)* sum(1/StdErr2[year,]^2, na.rm=TRUE)))    
      } else {
        sumse2 = sum(StdErr2[year,1:Nsurveys]^-2,na.rm=TRUE) 
        WAPW[year,] = c((StdErr2[year,1:Nsurveys]^-2)/sumse2,0)
        .YOutput$LogPred[year] = sum(WAPW[year,]*Pred[year,],na.rm=TRUE)      
        sumdiff = sum((Pred[year,1:Nsurveys]-.YOutput$LogPred[year])^2/(StdErr2[year,1:Nsurveys]^2),na.rm=TRUE)         
        tmp = Prediction[year,-(1:2)]
        nsury = length(tmp[!is.na(tmp)]) #Number of surveys positively reporting in this year
        sinvsq = sum((StdErr2[year,1:Nsurveys]^-2),na.rm=TRUE) 
        .YOutput$ISE[year] = (sinvsq)^-0.5
        .YOutput$ESE[year] = sqrt(sumdiff / ((nsury-1)* sum(1/StdErr2[year,1:Nsurveys]^2, na.rm=TRUE)))   
      }
     
    } else {
      warning(paste("xxxxxxx No survey data for",Prediction$Year[year]))
      .YOutput$LogPred[year] = NA
      .YOutput$ESE[year] = NA
      .YOutput$ISE[year] = NA
    }  
    .YOutput$Year[year] = Prediction$Year[year]
  }
  .YOutput$WAPred = exp(.YOutput$LogPred)-1
  .YOutput$VarRatio = (.YOutput$ESE/.YOutput$ISE)^2
  
#Access functions  
  list(
######  
    plot = function() 
    {
      predy = Index*Gradient+Intercept
      par(mfrow=c(Nprediction+1,Nsurveys+ifelse(shrink,3,2)),mar=c(0,0,0,0),oma=c(0,0,0,0)) #
      for (year in 1:Nprediction)
      {
        #yrange = c(0.8,1.2) * range(c(log(Calibration[,2]+1),predy[year,predy[year,]>0]),na.rm=TRUE) 
        yrange = c(0.8,1.2) * range(c(log(Calibration[,2]+1),predy[predy>0]),na.rm=TRUE) 
        yrange[2] = min(yrange[2],12*log(10))   ###greater than exp(700) numbers become infite
        yrange[1] = max(yrange[1],0.00000001)
        #browser()
        yy = (yrange[2]-yrange[1])*(0:1000)/1000 + yrange[1]
        if (!is.na(.YOutput$LogPred[year]))
        {
          plot(-dnorm(yy,.YOutput$LogPred[year],.YOutput$ISE[year]),exp(yy)-1,  
               axes=FALSE, ylim = exp(yrange)-1, frame.plot=TRUE,type='l',log='y')      
          yy = yy[yy<=.YOutput$ISE[year]+.YOutput$LogPred[year]&yy>=.YOutput$LogPred[year]-.YOutput$ISE[year]]
          polygon(c(-dnorm(yy,.YOutput$LogPred[year],.YOutput$ISE[year]),0,0),
                  exp(c(yy,.YOutput$LogPred[year]+.YOutput$ISE[year],.YOutput$LogPred[year]-.YOutput$ISE[year]))-1, 
                  col="green") 
          axis(4,las=1)
        } else {
          #No valid surveys for the year - can't do much!
          plot(1,1,axes=FALSE,type="n",frame.plot=TRUE,xlab="",ylab="")
          text(1,1,"No valid\nsurveys",col="red")      
        }
        plot(1,1,axes=FALSE,type="n",xlab="",ylab="")
        text(1,1,Prediction$Year[year])

        for (s in 1:Nsurveys)
        {
          xrange = c(0.8,1.2) * range(c(Calibration[,s+2],Prediction[,s+2]),na.rm=TRUE)
          if(!is.na(Gradient[year,s]))
          {
            
            #Plot points
            plot(Calibration[,s+2],Calibration$VPA,  
                  yaxt= ifelse(s==1,"s","n"), xaxt= ifelse(year==Nprediction,"s","n"),
                  xlim = xrange, ylim = exp(yrange)-1, log='y', las=1, frame.plot=TRUE)
            #Plot regression lines            

            xse = (0:100)*(xrange[2]-xrange[1])/100+xrange[1]
            stds = StdErr[year,s]*sqrt(sumw[year,s]/(sumw[year,s]-2)) *
                     sqrt(1+(1/sumw[year,s])+ (xse - xbar[year,s])^2/.sumsqx[year,s])  
            stdx = StdErr[year,s]*sqrt(sumw[year,s]/(sumw[year,s]-2)) *
                     sqrt(1+(1/sumw[year,s])+ (Index[year,s] - xbar[year,s])^2/.sumsqx[year,s])  
            yse = xse*Gradient[year,s]+Intercept[year,s]
            points(xse,exp(yse       )-1,type='l',col='red')            
            points(xse,exp(yse + stds)-1,type='l',col='red',lty="dotted")
            points(xse,exp(yse - stds)-1,type='l',col='red',lty="dotted")            
            #Plot prediction
            x = Index[year,s]
            y = Index[year,s]*Gradient[year,s]+Intercept[year,s]
            points(c(x,-999),exp(c(y,y))-1,type='l',col='green')
            points(c(x,x),exp(c(0.001,y+stdx))-1,type='l',col='green')
            points(c(x,-999),exp(c(y+stdx,y+stdx))-1,type='l',col='green',lty="dotted")
            points(c(x,-999),exp(c(y-stdx,y-stdx))-1,type='l',col='green',lty="dotted")            
          } else {
            if (is.na(Prediction[year,s+2]))
            {
              plot(Calibration[,s+2],Calibration$VPA,  
                    yaxt= ifelse(s==1,"s","n"), xaxt= ifelse(year==Nprediction,"s","n"),
                    xlim = xrange, ylim = exp(yrange)-1, log='y', las=1, frame.plot=TRUE)
                    text(mean(xrange),exp(mean(yrange))-1,"No survey\nthis year",col='red') 
                    #What if there is a survey, but rejected for not having enough points?
            } else {
              plot(Calibration[,s+2],Calibration$VPA,  
                    yaxt= ifelse(s==1,"s","n"), xaxt= ifelse(year==Nprediction,"s","n"),
                    xlim = xrange, ylim = exp(yrange)-1, log='y', las=1, frame.plot=TRUE)
                    text(mean(xrange),exp(mean(yrange))-1,"Not enough\npoints for\nprediction",col='red') 
                    #What if there is a survey, but rejected for not having enough points?            
            }
            
          }
        }
        if (shrink)
        {
          plot(rep(0,Ncalibration),Calibration$VPA,  
             yaxt= "n", xaxt= "n", ylim = exp(yrange)-1, log='y', las=1, frame.plot=TRUE)
          abline(exp(Intercept[year,Nsurveys+1])-1,0,col='red',untf=TRUE)
          abline(exp(Intercept[year,Nsurveys+1]+StdErr2[year,Nsurveys+1])-1,0,col='red',lty="dotted",untf=TRUE)
          abline(exp(Intercept[year,Nsurveys+1]-StdErr2[year,Nsurveys+1])-1,0,col='red',lty="dotted",untf=TRUE)
        }
      }
      plot(1,1,axes=FALSE,type="n",xlab="",ylab="")
      text(1,1,"ISE")
      plot(1,1,axes=FALSE,type="n",xlab="",ylab="")
      for (s in 1:Nsurveys) 
      {  
        plot(1,1,axes=FALSE,type="n",xlab="",ylab="")
        text(1,1,paste("log(",colnames(Prediction)[s+2],")",sep=""))

      }
    },
######    
    input = function() .input,
######    
    output =  function()
    {
      list(Gradient=Gradient,Intercept=Intercept,StdErr=StdErr,Rsquare=Rsquare,
           Npoints=Npoints,Index=Index,Pred=Pred,StdErrPred=StdErr2,WAPW=WAPW,Years=.YOutput)
    },
######    
    summary = function() .YOutput,
######    
    printout = function(filename)
    {
      opfile = file(filename,"w")
      cat(OPhead, '\n\n', file = opfile)
      for (y in 1:Nprediction)
      {
        cat("yearclass = ", Prediction$Year[y], "\n\n", file = opfile,append = TRUE,sep="")
        cat("          I-----------Regression----------I  I-----------Prediction---------I", "\n\n", file = opfile,append = TRUE,sep="")        
        cat(" Survey/  Slope  Inter-   Std  Rsquare  No.  Index Predicted   Std     WAP\n", file = opfile,append = TRUE,sep="") 
        cat(" Series           cept   Error          Pts  Value   Value    Error   Weights", "\n\n", file = opfile,append = TRUE,sep="")     
        for (s in 1:Nsurveys)
        {
        
          cat(sprintf2(8, " %7s",colnames(Prediction)[s+2]), file = opfile,append = TRUE,sep="") 
          cat(" ",sprintf2(6,"%6.2f",Gradient[y,s]), file = opfile,append = TRUE,sep="")                   
          cat(" ",sprintf2(6,"%6.2f",Intercept[y,s]), file = opfile,append = TRUE,sep="")                   
          cat(" ",sprintf2(6,"%6.2f",StdErr[y,s]), file = opfile,append = TRUE,sep="")                                       
          cat(" ",sprintf2(6,"%6.3f",Rsquare[y,s]), file = opfile,append = TRUE,sep="")      
          cat(" ",sprintf2(6,"%6i",as.integer(Npoints[y,s])), file = opfile,append = TRUE,sep="")   
          cat(" ",sprintf2(6,"%6.2f",Index[y,s]), file = opfile,append = TRUE,sep="")                             
          cat(" ",sprintf2(7,"%7.2f",Pred[y,s]), file = opfile,append = TRUE,sep="")                
          cat(" ",sprintf2(8,"%8.3f",StdErr2[y,s]), file = opfile,append = TRUE,sep="")                
          cat(" ",sprintf2(8,"%8.3f",WAPW[y,s]),"\n", file = opfile,append = TRUE,sep="")                                                                    
        }            
        cat("\n                                        VPA Mean =", file = opfile,append = TRUE,sep="")   
        cat(" ",sprintf2(7,"%7.2f",Pred[y,Nsurveys+1]), file = opfile,append = TRUE,sep="")                
        cat(" ",sprintf2(8,"%8.3f",StdErr2[y,Nsurveys+1]), file = opfile,append = TRUE,sep="")                
        cat(" ",sprintf2(8,"%8.3f",WAPW[y,Nsurveys+1]),"\n\n\n", file = opfile,append = TRUE,sep="")               
      }
      cat("\n\n\n\n Year     Weighted      Log     Int     Ext     Var     VPA      Log\n", file = opfile,append = TRUE,sep="")
      cat(" Class     Average      WAP     Std     Std    Ratio             VPA\n", file = opfile,append = TRUE,sep="")
      cat("          Prediction           Error   Error\n\n", file = opfile,append = TRUE,sep="")
      for (y in 1:Nprediction)
      {
        cat(" ",sprintf2(4,"%4i",.YOutput[y,1]), file = opfile,append = TRUE,sep="")
        cat(" ",sprintf2(11,"%11i",round(.YOutput[y,2])), file = opfile,append = TRUE,sep="")
        cat(" ",sprintf2(9,"%9.2f",.YOutput[y,3]), file = opfile,append = TRUE,sep="")        
        cat(" ",sprintf2(7,"%7.2f",.YOutput[y,4]), file = opfile,append = TRUE,sep="")
        cat(" ",sprintf2(7,"%7.2f",.YOutput[y,5]), file = opfile,append = TRUE,sep="")
        cat(" ",sprintf2(8,"%8.2f",.YOutput[y,6]),'\n', file = opfile,append = TRUE,sep="")
      }
      close(opfile)     
    }
  )
}


