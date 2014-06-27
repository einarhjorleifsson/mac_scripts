# Things to try:
#  rp.grid to improve layout
#  rp.doublebutton for taper and tapert, minpoints
#  rp.slider for minse
#  array of rp.checkbox for prelogged data
#  removing objects if file read doesn't work, or closing and reloading window

RCTinterface = function()
{

library(rpanel)
source(".\\RCT3v4a.r")

RCTchoosefile = function(panel)
{
  panel$Parameters$filename = choose.files(panel$Parameters$filename,caption = "Select data file",multi=FALSE)

  if (length(panel$Parameters$filename)==0) 
  { 
    getfile()
  } else {
    if (panel$Parameters$filename[1]=="") 
    {
      getfile()
    } else {
      #mainpanel <- rp.control(title = "RCT Interface", size = c(400,900), aschar=FALSE, RCT=list(), Parameters=panel$Parameters)
      mainpanel <- rp.control(title = "RCT Interface", aschar=FALSE, RCT=list(), Parameters=panel$Parameters)
      RCTsetup(mainpanel)
    }
  }

  panel
}

RCTsetup = function(panel)
{
  rp.radiogroup(panel, Parameters$shrink,values=c("y","n"),c("Shrinkage to mean","No Shrinkage"), title="Shrinkage", action = RCTplot)
  rp.radiogroup(panel, Parameters$regression,values=c("c","p"),c("Calibration","Prediction"), title="Regression type", action = RCTplot)
  rp.doublebutton(panel, Parameters$taper, step=1, "Power of taper", action=RCTplot,initval=3, range=c(0,5), showvalue=TRUE)
  rp.doublebutton(panel, Parameters$tapert, step=1, "Taper time period", action=RCTplot,initval=10, range=c(1,NA), showvalue=TRUE)
  rp.doublebutton(panel, Parameters$minpoints, step=1, "Minimum Points", action=RCTplot,initval=3, range=c(3,NA), showvalue=TRUE)
  rp.slider(panel, Parameters$minse, 0, 1, action = RCTplot,"Minimum SE",initval=0.2, showvalue=TRUE)
  rp.button(panel, RCTprintout, "Write output to file")
  rp.button(panel, RCToutput, "Full output to terminal")
  rp.button(panel, RCTsummary, "Summary output to terminal")
  rp.button(panel, RCTinput, "Input data to terminal")
  rp.button(panel, getfile, "Different Data", quitbutton=TRUE)
  rp.button(panel, noAction, "Exit", quitbutton=TRUE)
#  rp.radiogroup(panel, Parameters$shrink,values=c("y","n"),c("Shrinkage to mean","No Shrinkage"), title="Shrinkage", action = RCTplot,pos=c(0,0,200,100))
#  rp.radiogroup(panel, Parameters$regression,values=c("c","p"),c("Calibration","Prediction"), title="Regression type", action = RCTplot,pos=c(200,0,200,100))
#  rp.doublebutton(panel, Parameters$taper, step=1, "Power of taper", action=RCTplot,initval=3, range=c(0,5), showvalue=TRUE,pos=c(0,100,200,100))
#  rp.doublebutton(panel, Parameters$tapert, step=1, "Taper time period", action=RCTplot,initval=10, range=c(1,NA), showvalue=TRUE,pos=c(200,100,200,100))
#  rp.doublebutton(panel, Parameters$minpoints, step=1, "Minimum Points", action=RCTplot,initval=3, range=c(3,NA), showvalue=TRUE,pos=c(300,0,200,100))
#  rp.slider(panel, Parameters$minse, 0, 1, RCTplot,"Minimum SE",initval=0.2, showvalue=TRUE,pos=c(300,200,200,100))
#  rp.button(panel, RCTprintout, "Write output to file",pos=c(400,50,200,80))
#  rp.button(panel, RCToutput, "Full output to terminal",pos=c(500,50,200,80))
#  rp.button(panel, RCTsummary, "Summary output to terminal",pos=c(600,50,200,80))
#  rp.button(panel, RCTinput, "Input data to terminal",pos=c(700,50,200,100))
#  rp.button(panel, getfile, "Different Data", quitbutton=TRUE,pos=c(800,50,200,80))
  panel = RCTplot(panel)
  panel
}

#test = function(panel)
#{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
#  panel2 <- rp.control("RCT Interface", aschar=FALSE, firsttime=TRUE, RCT=list(), Parameters=Parameters)
#  panel
#}

RCTprintout = function(panel)
{

  opfile = choose.files("c:\\rct3\\data\\output.txt")
  if (length(opfile)>0) panel$RCT$printout(opfile)
  panel
}

RCTsummary = function(panel)
{
  print(panel$RCT$summary())
  panel
}

RCTinput = function(panel)
{
  print(panel$RCT$input())
  panel
}

RCToutput = function(panel)
{
  print(panel$RCT$output())
  panel
}

RCTplot = function(panel)
{
  panel$RCT = RCT3(filename   = panel$Parameters$filename, 
                   old        = panel$Parameters$old, 
                   logged     = panel$Parameters$logged, 
                   minpoints  = panel$Parameters$minpoints, 
                   regression = panel$Parameters$regression, 
                   shrink     = (panel$Parameters$shrink=="y"),
                   taper      = panel$Parameters$taper, 
                   tapert     = panel$Parameters$tapert, 
                   minse      = panel$Parameters$minse, 
                   configfile = panel$Parameters$configfile)
  if (is.null(panel$RCT))
  {
    panel$RCT = RCT3(filename   = panel$Parameters$filename, 
                     old        = !panel$Parameters$old, 
                     logged     = panel$Parameters$logged, 
                     minpoints  = panel$Parameters$minpoints, 
                     regression = panel$Parameters$regression, 
                     shrink     = (panel$Parameters$shrink=="y"),
                     taper      = panel$Parameters$taper, 
                     tapert     = panel$Parameters$tapert, 
                     minse      = panel$Parameters$minse, 
                     configfile = panel$Parameters$configfile)
  if (is.null(panel$RCT))
  {
    rp.messagebox("File type not recognised")
    return(panel)
  }
  
  }
  if (is.null(dev.list())) x11(14,10)
  panel$RCT$plot()
  panel 
}

noAction = function(panel)
{
  panel
}


                  


#
#ple7arec = RCT3("c:\\rct3\\data\\ple7arec.csv")
#ple7arec$input()
#ple7arec$plot()
#ple7arec$summary()
#ple7arec$output()
#ple7arec$printout("c:\\rct3\\data\\ple7arecOutput.txt")
#should match  ple7arc3.txt

getfile = function(panel=NA)
{
  graphics.off()
  Parameters = list(filename = "c:\\rct3\\data\\ple7arec.csv",
                     old = TRUE, 
                     logged = NA, 
                     minpoints = 3, 
                     regression='c', 
                     shrink="y",
                     taper=3, 
                     tapert=20, 
                     minse = 0, 
                     configfile = NA)
  filepanel <- rp.control("RCT Interface", size=c(200,100),aschar=FALSE, Parameters=Parameters)
  rp.button(filepanel, RCTchoosefile, title = "Choose input file", quitbutton=TRUE,pos = c(25,25,150,50))
  panel
}

getfile()
  5

}