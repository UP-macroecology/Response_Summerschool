# Example RangeShiftR workflow
# RTG Response summer school 2021
# (c) Damaris Zurell, Anne Malchow, Univ. Potsdam
# More elaborate materials: https://rangeshifter.github.io/RangeshiftR-tutorials/ and https://damariszurell.github.io/EEC-QCB/

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
if(!file.exists("models/GettingStarted")) {
  dir.create("models/GettingStarted", showWarnings = TRUE) }


# The standard workflow of RangeShiftR is to load input maps from ASCII raster files and to write all simulation output into text files. Therefore, the specified working directory needs to have a certain folder structure: It should contain 3 subfolders named "Inputs", "Outputs" and "Output_Maps".

# relative path from working directory:
dirpath = "models/GettingStarted/"

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

# RangeShiftR can either import a map from an ASCII raster file in the ‘Inputs’ folder or generate a random map to use in the simulation.

# For first illustration, we use an artifical landscape
land <- ArtificialLandscape(Resolution = 10,  # in meters
                            K_or_DensDep = 1500,  # ~ 15 inds/cell
                            propSuit = 0.2,
                            dimX = 129, dimY = 257, 
                            fractal = T, hurst = 0.3,
                            continuous = F)

#---------------
# DEMOGRAPHY
#---------------

# The Demography module contains all the local population dynamics of your simulated species. Generally there are two types, Unstructured model / non-overlapping generations, and Stage-structured model / overlapping generations

# This is how we can define an unstructured model:
demo <- Demography(Rmax = 2.2, ReproductionType = 1, PropMales = 0.45)

# Alternatively, a stage-structure model might look like this:
stg <- StageStructure(Stages = 3,
                      TransMatrix = matrix(c(0,1,0,5.7,.5,.4,3.4,0,.9),nrow = 3),
                      FecDensDep = T,
                      SurvDensDep = T)
demo <- Demography(StageStruct = stg, ReproductionType = 1, PropMales = 0.45)

# RangeShiftR provides a number of useful functions to explore the model set-up. For example, we can plot the rates from the transition matrix:
plotProbs(stg)


#---------------
# DISPERSAL
#---------------

# The dispersal process is modelled wih three sub-processes (see the schematic figure above): Emigration(), Transfer() and Settlement().

disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.2), 
                   Transfer   = DispersalKernel(Distances = 50),
                   Settlement = Settlement() )

# Again, we can use the plotProbs() function to plot the functional relationships, for example for the transfer phase
plotProbs(DispersalKernel(Distances = 50))


#---------------
# INITIALISATION
#---------------

# In order to control the initial distribution of individuals in the landscape at year 0, we set initialisation rules. We choose to initialise 3 individuals per habitat cell. Additionally, since we define a stage-structured model, we have to specify the initial proportion of stages:

init <- Initialise(FreeType = 0, 
                   NrCells = 2250,
                   InitDens = 2, 
                   IndsHaCell = 3, 
                   PropStages = c(0,0.7,0.3))
init


#-----------------
# SIMULATION
#-----------------

# This module is used to set general simulation parameters (e.g. simulation ID, number of replicates, and number of years to simulate) and to control output types (plus some more specific settings).

sim <- Simulation(Simulation = 2,
                  Years = 50,
                  Replicates = 2,
                  OutIntPop = 50)

#-----------------
# PARAMETER MASTER
#-----------------

# After all settings have been made in their respective modules, we are ready to combine them to a parameter master object, which is needed to run the simulation.

s <- RSsim(simul = sim, land = land, demog = demo, dispersal = disp, init = init)

# Alternative notation:
# s <- RSsim() + land + demo + disp + sim + init

# We can check the parameter master (or any single module) for potential parameter conflicts:
validateRSparams(s)


#-----------------
# RUN SIMULATION
#-----------------

# Once the parameter master has been defined, we can run the simulations in the specified RS directory.
RunRS(s, dirpath)


#-----------------
# PLOT RESULTS
#-----------------

# All results are stored in the Outputs folder. 
# RangeShiftR provides some in-built functions to access and plot these results. 
# First, we plot the abundance and occupancy time series:
range_df <- readRange(s, dirpath)

# plot mean abundance and mean occupancy along with different replicate runs (n=2):
par(mfrow=c(1,2))
plotAbundance(range_df)
plotOccupancy(range_df)

# plot mean abundance and mean occupancy along with standard deviation across replicates (n=2):
par(mfrow=c(1,2))
plotAbundance(range_df, sd=T, replicates = F)
plotOccupancy(range_df, sd=T, replicates = F)


# Based on the file outputs, you can also post-process the data yourself and make more complex output plots.
# for example, plot the spatial distribution of abundance:

# read population output file into a dataframe
pop_df <- readPop(s, dirpath)

# Not all years have the same number of populated and thus listed cells. For stacking, we set a common extent with the values used in the landscape module:
ext <- c(0,1290,0,2570)
res <- 10

# Make stack of different raster layers for each year and for only one repetition (Rep==0):
pop_wide_rep0 <- reshape(subset(pop_df,Rep==0)[,c('Year','x','y','NInd')], timevar='Year', v.names=c('NInd'), idvar=c('x','y'), direction='wide')

# use raster package to make a raster from the data frame
stack_years_rep0 <- rasterFromXYZ(pop_wide_rep0)
names(stack_years_rep0) <- c('Year.0', 'Year.50')
spplot(stack_years_rep0, zlim = c(0,10))

