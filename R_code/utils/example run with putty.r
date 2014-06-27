##########################################################################################################
# 
#  North Sea Flatfish Management Plan:  Evaluation of proposed changes
#  D. Miller, A. Coers, J.J. Poos
#  For AGVAMPS
# 
#  Last updated: 18/09/2012
# 
##########################################################################################################
rm(list=ls(all=TRUE))
memory.limit(size=4000)
### ------------------------------------------------------------------------------------------------------
###   1. Load libraries, set path, source functions
### ------------------------------------------------------------------------------------------------------
### set path:
putty <- T 
if (!putty) myPath <- "N:/Projecten/Flatfish strategy evaluation 2012/" else myPath <- "/media/n/Projecten/Flatfish\ strategy\ evaluation\ 2012/"
# Requires R 2.13.2 with appropriate FLR versions
# Install packages  (FLCore, FLAssess, FLXSA)
# install.packages(repos="http://flr-project.org/R")
# Load libraries
library(FLCore); library(FLAssess); library(FLXSA)
#library(FLEDA); library(FLSTF); 
library(splines); library(mvtnorm);library(MASS); library(boot);

### ------------------------------------------------------------------------------------------------------
###   0. Load Scenarios, set up version saving
### ------------------------------------------------------------------------------------------------------
# Run name (scenario Set)
scenSet <- "1_Comparison"
# Create results directory
setwd(paste(myPath,"3_Data/3-2_Finale_data/", sep=""))
if (!putty) shell(paste("mkdir",scenSet))

# Choose scenarios to run (see "FF_Scenarios_(scenSet).csv"): 
scenarios <- c(1:5)       
# Scenario details   
potScens <- read.csv(paste(myPath, "8_Code/Scenarios/FF_Scenarios_",scenSet,".csv",sep=""), header=T)
verList <- list()


### ------------------------------------------------------------------------------------------------------
###   1. Load data and get SR relationships
### ------------------------------------------------------------------------------------------------------
### data settings:
### Years
XSAassYear <- 2012     # latest XSA assessment year
#SCAassYear <- 2012     # latest SCA assessment year
#scaStart <- 1985    # first year of SCA output to use

# Dimensions
runs            <-  200                           # number of iterations (100)
startYear       <-  2003                        # 1st year of simulation   (= first year with observed effort data)
finalYear       <-  2026                          # last year of simulation  (i.e. have population at the start)

### Stocks
# NOTE: need Lowestoft format data files in the directoy(s): paste(myPath, "3 Data/3-1\ ruwe\ data/", stockName, "/,assyear,/", sep="")
ffSpp           <- c("PLE", "SOL")                    # Flatfish species
startPts        <- c("XSA","best","worst")#,"SCA","SCA_5","SCA_95")    # Starting points: XSA, best case, worst case , SCA percentiles
bestSP <- c(0.25,0.05)                           # Mean and SD for log normal distribution around multiplication factor for best case
worstSP <- c(-0.25,0.05)                         # Mean and SD for log normal distribution around multiplication factor for worst case
numCohorts <- 6                                  # number of latest cohorts to adjust

# Data been loaded previously?  (this is TRUE every time after the first time, unless changes are made to '01 read data.r')
# Re-run Data if any of the above has been changed
# MUST BE SET TO TRUE WHEN RUNNING ON PUTTY (i.e. load data running on computer first)
dataLoaded <- T
# Run data once:
if (!dataLoaded) source(paste(myPath,"8_Code/01_Read_data_2012.r", sep=""))
# Change myPath (will have been saved differently in run data (could remove there)
putty <- T
if (!putty) myPath <- "N:/Projecten/Flatfish strategy evaluation 2012/" else myPath <- "/media/n/Projecten/Flatfish\ strategy\ evaluation\ 2012/"

# Check Data Load:
#indRes[["XSA"]][["SOL"]][[1]]
  
### ------------------------------------------------------------------------------------------------------
###   2. Start Scenario loop
### ------------------------------------------------------------------------------------------------------
for (Sc in 1:length(scenarios)) {
#Sc <- 1
#Sc <- 2

# Remove past data (except a few objects)
blah <- ls()
keepObjs <- c("Sc","myPath","scenarios","potScens","verList","scenSet","obs_error")
for (ko in 1:length(keepObjs)) blah <- blah[-(which(blah==keepObjs[ko]))]
rm(list=blah); rm(blah)

# Save verList
save.image(file = paste(myPath,"3_Data/3-2_Finale_data/",scenSet,"/verList_",scenSet,".RData", sep=""))
                  
### ------------------------------------------------------------------------------------------------------
###   2a. Load data and get SR relationships
### ------------------------------------------------------------------------------------------------------
load(paste(myPath, "3_Data/3-1_Ruwe_data/inputData_",scenSet,".RData",sep=""))
putty <- T
if (!putty) myPath <- "N:/Projecten/Flatfish strategy evaluation 2012/" else myPath <- "/media/n/Projecten/Flatfish\ strategy\ evaluation\ 2012/"

### ------------------------------------------------------------------------------------------------------
###   2b. Set model options
### ------------------------------------------------------------------------------------------------------
### Simulation settings:
# Stock
startPt <- ac(potScens[scenarios[Sc],"startPt"])     # "XSA", "SCA", "SCA_5" or "SCA_95"
SRtype <- list()                                     # SR relationship(s) - specify for each of ffSpp: geomean, bevholt, minrec (minimum of timeseries)   #NO RICKER!
SRtype[["PLE"]]   <- ac(potScens[scenarios[Sc],"SRtype"])                     
SRtype[["SOL"]]   <- ac(potScens[scenarios[Sc],"SRtype"])  
numYmn          <- 5                                 # Number of years to average/resample for future values e.g. selectivity, weights etc.

# Fleets
fltList       <- c("NL_BT2","NL_Other","Other")     # Fleets: choose from c("NL_BT2","NL_Other","Other")  - # Add option for total
relStabNL <- list()
relStabNL[["PLE"]]   <- "mean"  #0.37 or "mean"      # NL Fleet relative stability for PLE ('Other' Fleet rel. stab. = 1 - this)   (values from EU TAC allocations - or used past observed landings ratios (0.4)?)
relStabNL[["SOL"]]   <- "mean"  #0.75 or "mean"      # NL Fleet relative stability for SOL ('Other' Fleet rel. stab. = 1 - this)   (values from EU TAC allocations - or used past observed landings ratios (0.7)?)
stkEffort <- ac(potScens[scenarios[Sc],"stkEffort"]) # Which stock determines effort limit: "most" (maximum allowed species effort), "both" (no over quota, species specific allowable effort) or "least" (minimum allowed species effort)
upperLim  <- ac(potScens[scenarios[Sc],"upperLim"])  # Catch with all Effort ("effort") or stop once TAC is caught ("TAC")
impEffort <- ac(potScens[scenarios[Sc],"impEffort"]) # Implementation effort(over quota): "most" (maximum allowed species effort), "both" (no over quota, species specific allowable effort) or "least" (minimum allowed species effort)
pastQ     <- ac(potScens[scenarios[Sc],"pastQ"])     # "obs" or "gen"              # Q for past years i.e. 2003-2008 (as observed each year or generated from mean and SD over the whole period)
techCreep <- list()                                  # technical creep values (as decimals) - specify for each of ffSpp
techCreep[["PLE"]]   <- as.numeric(potScens[scenarios[Sc],"PLE_TC"])    # Rijnsdorp et al 2006: 2.8% for sole, 1.6% for plaice
techCreep[["SOL"]]   <- as.numeric(potScens[scenarios[Sc],"SOL_TC"])  
#ADD Future q options here if needed
FvEff_rel <- list()
FvEff_rel[["Lan"]] <- 1                             # Relationship between landings q and F (1=linear)
FvEff_rel[["Dis"]] <- 1                             # Relationship between discards q and F (1=linear)

## check to see if fleet selection covers total effort
if (sum(apply(fltEff[fltList,],2,sum,na.rm=T))<sum(fltEff["Total",])) print("Fleet selection does not include ALL effort!")
if (sum(apply(fltEff[fltList,],2,sum,na.rm=T))>sum(fltEff["Total",])) print("Fleet selection exceeds ALL effort! i.e. there is overlap")

# Harvest control rule values
HCR_args <- list(name = "676/2007",                   
                 FsqMethod = as.character(potScens[scenarios[Sc],"FsqMethod"]),   # ICES = prev. year rescaled; doubICES = prev. year rescaled, less reduction; TAC = use TAC.
                 maxTACchange = as.numeric(potScens[scenarios[Sc],"HCR_TACchng"]),
                 article18 = T,    # Allows for greater than 15% reduction (guessed at 25% or % below) if stock outside safe bio limits.
                 Fred = as.numeric(potScens[scenarios[Sc],"HCR_Fred"]),
                 effRed = "rel",    #can be "fixed" (always 10%) or "rel" reduce relative to amount above Ftar (max 10%)
                 effCap  =  potScens[scenarios[Sc],"effCap"],  
                 effCapLev = 28307876, 
                 maxTAEchange = as.numeric(potScens[scenarios[Sc],"HCR_TAEchng"]),
                 Bpa     = list(PLE=230000,SOL=35000),
                 Blim    = list(PLE=160000,SOL=25000),
                 Fpa     = list(PLE=0.6,SOL=0.4),
                 Flim    = list(PLE=0.74,SOL=NA),
                 Ftar    = list(PLE=as.numeric(potScens[scenarios[Sc],"HCR_FtarPLE"]), SOL=as.numeric(potScens[scenarios[Sc],"HCR_FtarSOL"])),
                 Fages  = list(PLE=2:6, SOL=2:6 )
                 )                                  

                 
### ------------------------------------------------------------------------------------------------------
###   3. Run simulation
### ------------------------------------------------------------------------------------------------------
# NOTE: files "Flatfish TACs.csv" (historic TACs for each of "ffSpp") and "fleet_effort.csv" (historic effort for each of "natFleets") in directory: paste(myPath,"3 Data/3-1 Ruwe data/",sep="")
#set.seed(0)
# Record unique version number
#VerNum          <- paste(ac(potScens[scenarios[Sc],"Number"]),"-",ac(potScens[scenarios[Sc],"Description"]),"_",startPt,"_",SRtype[1],"_",if(HCR_args$Ftar[[1]]!=0.2) "altFtar_",if(HCR_args$maxTACchange!=0.15) paste("TACchng-",HCR_args$maxTACchange,"_",sep=""),if(HCR_args$Fred!=0.1) paste("Fred-",HCR_args$Fred,"_",sep=""),runs,"reps-",format(Sys.time(), "%d%b%Y_%Hh%M"), sep="")
VerNum          <- paste(ac(potScens[scenarios[Sc],"Number"]),"-",ac(potScens[scenarios[Sc],"Description"]),"_",startPt,"_",SRtype[1],"_",if(HCR_args$Ftar[[1]]!=0.2) "altFtar_",if(HCR_args$effCap) "effCap_",runs,"reps-",format(Sys.time(), "%d%b%Y_%Hh%M"), sep="")

# Run model
source(paste(myPath,"8_Code/02_Model_2012.r", sep=""))

### ------------------------------------------------------------------------------------------------------
###   4. Save scenario
### ------------------------------------------------------------------------------------------------------
save.image(file = paste(myPath,"3_Data/3-2_Finale_data/",scenSet,"/",VerNum,".RData", sep=""))

# Save VerNum in list
verList[[ac(potScens[scenarios[Sc],"Description"])]] <- VerNum

# Timer
print(paste("End of scenario: ",Sc,"/",length(scenarios),", Estimated time left: ",round((End-SttSim)["elapsed"]/60,2)*(length(scenarios)-Sc),"min",sep=""))

  } # end of Scenario loop

# Save last verList
# Remove past data (except a few objects)
blah <- ls()
keepObjs <- c("Sc","myPath","scenarios","potScens","verList","scenSet","obs_error")
for (ko in 1:length(keepObjs)) blah <- blah[-(which(blah==keepObjs[ko]))]
rm(list=blah); rm(blah)
save.image(file = paste(myPath,"3_Data/3-2_Finale_data/",scenSet,"/verList_",scenSet,".RData", sep=""))
 
