# Instructions
# New screen
# Use '-S' and add a name afterwards
# e.g. screen -S (name)

# to exit screen (go back to main):
#ctrl-A + d (to exit the screen), i.e. Ctrl+A, release, then D

# to go back to the screen:
# screen -r (name)

# To list available screens:
# screen -ls

# To reattach screen
#  screen -DR (name)

# To kill screen
# kill (name)

# To see computer performance:
# htop
# close htop with F10


# For putty
## Run Main
#source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012.r")
## Run stats
#source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/04b_Stats_for_MSE_2012.r")
## Run plots
#source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/04_Plots_for_MSE_2012.r")
#
  
## Base case scenarios
screen -S TEST
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_00_TEST.r")
  
## Base case scenarios
screen -S 1_Comparison
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_1_Comparison.r")

## Base case scenarios short
screen -S Quickie
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_1_Comparison_short.r")

## Effcap scenarios
screen -S 2_BestWorst
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_2_BestWorst.r")

## Effcap scenarios
screen -S 3_Effort
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_3_Effort.r")

## New scenarios
screen -S 4_New
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_4_New.r")

## Implementation Effort scenarios
screen -S 2_ImpEffort
# Run Main
source("/media/n/Projecten/BO\ Flatfish\ strategy\ evaluation\ 2012/8_Code/00_Main_loop_2012_2_ImpEffort.r")
