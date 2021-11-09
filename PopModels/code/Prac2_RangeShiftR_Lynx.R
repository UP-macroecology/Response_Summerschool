# Example RangeShiftR: Lynx reintroduction
# RTG Response summer school 2021
# (c) Damaris Zurell, Anne Malchow, Univ. Potsdam
# Full tutorial available here: https://damariszurell.github.io/EEC-QCB/RS3_lynx.html
# More elaborate materials: https://rangeshifter.github.io/RangeshiftR-tutorials/ and https://damariszurell.github.io/EEC-QCB/


# This practical illustrates how spatially-explicit individual-based models like RangeShiftR can aid effective decision making for species reintroductions. As example, we re-implement a case study on the reintroduction of Eurasian lynx (Lynx lynx) to Scotland (Ovenden et al. 2019; https://doi.org/10.1016/j.biocon.2019.03.035). The original study was also based on RangeShifter, which makes reimplementation straight forward. We use the same parameters by and large but on a slightly coarser resolution. In line with Ovenden et al. (2019), we simulate lynx range expansion and population viability from different potential reintroduction sites.


#--------------------------------------------------------------------------------
#
#         Load packages
# 
#--------------------------------------------------------------------------------

library(RangeShiftR)  # RangeShiftR package for spatially-explicit eco-evolutionary modelling

library(raster)   # Manipulate geographic data


#--------------------------------------------------------------------------------
#
#         SET WORKING DRECTORY
# 
#--------------------------------------------------------------------------------

# Set your working directory to the workshop folder, e.g 
setwd('Response_Summerschool/PopModels')



#--------------------------------------------------------------------------------
#
#         PREPARE FOLDER STRUCTURE
# 
#--------------------------------------------------------------------------------

# We recommend a single models folder in your project folder and model subfolders for the different RangeShiftR practicals
# Set up folder "models" in working directory
if(!file.exists("models")) {
  dir.create("models", showWarnings = TRUE) }
# Set up subfolder "GettingStarted" within "models" folder
if(!file.exists("models/Lynx")) {
  dir.create("models/Lynx", showWarnings = TRUE) }


# The standard workflow of RangeShiftR is to load input maps from ASCII raster files and to write all simulation output into text files. Therefore, the specified working directory needs to have a certain folder structure: It should contain 3 subfolders named "Inputs", "Outputs" and "Output_Maps".

# relative path from working directory:
dirpath = "models/Lynx/"

# Create sub-folders (if not already existing)
if(!file.exists(paste0(dirpath,"Inputs"))) {
  dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Outputs"))) {
  dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Output_Maps"))) {
  dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE) }


#--------------------------------------------------------------------------------
#
#         RangeShiftR MODULES
# 
#--------------------------------------------------------------------------------


#---------------
# LANDSCAPE
#---------------



#---------------
# DEMOGRAPHY
#---------------




#---------------
# DISPERSAL
#---------------




#---------------
# INITIALISATION
#---------------




#-----------------
# SIMULATION
#-----------------



#-----------------
# PARAMETER MASTER
#-----------------




#-----------------
# RUN SIMULATION
#-----------------




#-----------------
# PLOT RESULTS
#-----------------



