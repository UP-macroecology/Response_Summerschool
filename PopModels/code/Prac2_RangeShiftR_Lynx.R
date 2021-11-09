# Example RangeShiftR: Lynx reintroduction
# RTG Response summer school 2021
# (c) Damaris Zurell, Anne Malchow, Univ. Potsdam
# Full tutorial available here: https://damariszurell.github.io/EEC-QCB/RS3_lynx.html
# More elaborate materials: https://rangeshifter.github.io/RangeshiftR-tutorials/ and https://damariszurell.github.io/EEC-QCB/
# Full RangeShifter manual: https://raw.githubusercontent.com/RangeShifter/RangeShifter-software-and-documentation/master/RangeShifter_v2.0_UserManual.pdf


# This practical illustrates how spatially-explicit individual-based models like RangeShiftR can aid effective decision making for species reintroductions. As example, we re-implement a case study on the reintroduction of Eurasian lynx (Lynx lynx) to Scotland (Ovenden et al. 2019; https://doi.org/10.1016/j.biocon.2019.03.035). The original study was also based on RangeShifter, which makes reimplementation straight forward. We use the same parameters by and large but on a slightly coarser resolution. In line with Ovenden et al. (2019), we simulate lynx range expansion and population viability from different potential reintroduction sites.


#--------------------------------------------------------------------------------
#
#         Load packages
# 
#--------------------------------------------------------------------------------

library(RangeShiftR)  # RangeShiftR package for spatially-explicit eco-evolutionary modelling

library(raster)   # Manipulate geographic data
library(rasterVis)    # Enhanced visualization of geographic data
library(ggplot2)    # Advanced visualisation
library(tidyverse)  # family of data science packages


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


# Download the files "LCM_Scotland_2015_1000.asc" and "woodland_patchIDs_1000.asc" from the Github repo and place them in the "models/Lynx/Inputs" folder.


#--------------------------------------------------------------------------------
#
#         RangeShiftR MODULES
# 
#--------------------------------------------------------------------------------


#---------------
# LANDSCAPE
#---------------

# The original lynx model was run on a 100 m landscape grid. Here, we use a spatial resolution of 1 km to speed up computations, for illustrative purposes.


#---------------
# Habitat classification:
#---------------


# Land cover map of scotland
landsc <- raster(paste0(dirpath, "Inputs/LCM_Scotland_2015_1000.asc"))
plot(landsc)


# Make prettier maps with proper legend:
# Plot land cover map and highlight cells with initial species distribution - option 2 with categorical legend:
landsc.f <- as.factor(landsc)

# add the land cover classes to the raster attribute table
(rat <- levels(landsc.f)[[1]])

rat[["land-use"]] <- c("salt water", "arable + horticulture", "freshwater", "built-up areas + gardens", "inland rock", "grasslands", "woodland", "supra-/littoral sediment", "marsh, swamp", "heath")
levels(landsc.f) <- rat

levelplot(landsc.f, margin=F, scales=list(draw=FALSE), col.regions=brewer.pal(n = 10, name = "Spectral"))


# We will run RangeShiftR as patch-based model and thus need a map of suitable patches (vs. unsuitable matrix). We have prepared a map of woodland patches from the above land cover map and will read this in. These woodland patches constitute the suitable breeding patches for the Eurasian lynx. In total, we have identified 39 suitable patches of varying size.

# Read in woodland patch files
patches <- raster(paste0(dirpath,'Inputs/woodland_patchIDs_1000.asc'))

# Plot the patches in different colours:
plot(patches,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 10, name = "Paired"),4)))

# In total, we have 39 patches of different size
values(patches) %>% table


#---------------
# Reintroduction patch:
#---------------

# Each woodland patch has a unique ID.

# Plot patches and add labels
plot(patches,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 10, name = "Paired"),4)))
patches_pol <- rasterToPolygons(patches, dissolve=T) # Makes spatial polygons
text(patches_pol,labels=patches_pol@data$woodland_patchIDs_1000)

# The following patch IDs roughly correspond to the reintroduction sites used in (Ovenden et al. 2019): 29 = Kintyre Peninsula; 35 = Kielder Forest; 15 = Aberdeenshire. Note that the patches are smaller and more isolated than in the original study.


#---------------
# Landscape parameters:
#---------------

# We can now define the landscape parameters. 
# For the patch-based model, we have to provide a (i) grid-based landscape map that contains the different habitat classes used to set the per-step mortality in the dispersal module, and (ii) the patch file used to define the suitable breeding patches. 
# Also, we have to provide the parameter K_or_DensDep for each habitat class in the landscape file. This parameter describes the strength of density dependence 1/b (Bocedi et al. 2020; Malchow et al. 2020). Ovenden et al. (2019) specified 1/b with 0.000285 ind/ha, meaning that an adult individual will defend a territory of c. 3509 ha or 35 km2.

land <- ImportedLandscape(LandscapeFile = "LCM_Scotland_2015_1000.asc",
                          PatchFile = "woodland_patchIDs_1000.asc", 
                          Resolution = 1000,
                          Nhabitats = 10,
                          K_or_DensDep = c(0, 0, 0, 0, 0, 0, 0.000285, 0, 0, 0)
)


#---------------
# DEMOGRAPHY
#---------------

# All demographic parameters are listed in Table 2 of Ovenden et al. (2019). The demographic model comprises three stages: first-year juveniles (0–12 months), non-breeding sub-adults (12–24 months) and breeding adults (>24months). To properly simulate natal dispersal, we have to set an additional dummy stage 0 (dispersing juveniles) in RangeShiftR. After dispersal, dispersing juveniles (stage 0) will be assigned to stage 1 (first-year juveniles). 

# Transition matrix
(trans_mat <- matrix(c(0,1,0,0,0,0, 0.53, 0,0, 0, 0, 0.63, 5,0, 0, 0.8), nrow = 4, byrow = F))  # stage 0: dispersing newborns; stage 1 : first-year juveniles; stage 2: non-breeding sub-adults; stage 3: breeding adults

# Define the stage structure
stg <- StageStructure(Stages = 4, # four stages including dispersing juveniles
                      TransMatrix = trans_mat,  # transition matrix
                      MaxAge = 17,  # maximum age
                      FecDensDep = T   # density dependence in fecundity
)

# Plot the vital rates against different density levels
plotProbs(stg)

# Define the demography module
demo <- Demography(StageStruct = stg, 
                   ReproductionType = 1, # simple sexual model
                   PropMales = 0.5) 



#---------------
# DISPERSAL
#---------------

# Ovenden et al. (2019) considered Eurasian lynx as poor dispersers with some sex bias in the dispersal probability. Specifically, the males were assigned a higher emigration probability compared to females.

emig <- Emigration(DensDep=T, StageDep=T, SexDep = T,
                   EmigProb = cbind(c(0,0,1,1,2,2,3,3),
                                    c(0,1,0,1,0,1,0,1),
                                    c(0.4, 0.9, 0, 0, 0, 0, 0, 0),
                                    c(10, 10, 0, 0, 0, 0, 0, 0), 
                                    c(1, 1, 0, 0, 0, 0, 0, 0)) ) # only emigration of juveniles, females higher than males


# The transfer phase of dispersal was modelled using a stochastic movement simulator (SMS). 
# Individuals move stepwise from cell to cell and the direction chosen at each step is determined by the land cover costs (specified for each habitat class), the species’ perceptual range (PR) and directional persistence (DP). 
transfer <- SMS(PR = 1, # Perceptual range in number of cells
                PRMethod = 2, # Harmonic mean used to quantify movement cost within perceptual range
                MemSize = 5, # number of steps remembered when applying directional persistence.
                DP = 5,  # directional persistence.
                Costs = c(100000, 30, 100, 1000, 1000, 10, 1, 10, 10, 7), # movement cost per habitat class
                StepMort = c(0.9999, 0.0002, 0.0005, 0.007, 0.00001, 0.00001, 0, 0.00001, 0.00001, 0.00001)) # per step mortality per habitat class


# The settlement parameters are identifcal for both sexes except that males have to find a female to settle. 
settle <- Settlement(StageDep = F,
                     SexDep = T,
                     Settle = cbind(c(0, 1), c(1.0, 1.0), c(-10, -10), c(1, 1)), # here no difference between sexes
                     FindMate = c(F, T), # only males need to find a female
                     DensDep = T,
                     MaxSteps = 500
)

# Define Dispersal module
disp <-  Dispersal(Emigration = emig,
                   Transfer = transfer,
                   Settlement = settle) 

# plot parameters
plotProbs(disp@Emigration)


#---------------
# INITIALISATION
#---------------

# Ovenden et al. (2019) implemented single-site reintroductions with 10 individuals per site, and multi-site reintroduction with 32 individuals across two sites.


#---------------
# Single-site reintroduction
#---------------

# We define one scenario where we reintroduce 10 individuals to the patch 29 on the Kintyre Pensinsula. First, we have to prepare a text file specifying the number of initial individuals per patch and per sex and stage.

# prepare dataframe for InitIndsFile
(init_df_29 <- data.frame(Year=0,Species=0,PatchID=29,Ninds=c(5,5),Sex=c(0,1),Age=3,Stage=3))

# We write the list of initial individuals into a file in the Inputs folder and then specify the initialisation module.

# write InitIndsFile to file
write.table(init_df_29, file=paste0(dirpath,'Inputs/InitInds_29.txt'), sep='\t', row.names=F, quote=F)

# Set initialisation
init_29 <- Initialise(InitType = 2,       # Initialise from initial individuals list file
                      InitIndsFile = 'InitInds_29.txt')


#---------------
# Multi-site reintroduction
#---------------

# In the multi-site reintroduction scenarios, Ovenden et al. (2019) reintroduced 18 lynx to Kintyre Peninsula and 14 lynx to Aberdeenshire.

# prepare dataframe for InitIndsFile
(init_df_29_15 <- data.frame(Year=0,Species=0,PatchID=rep(c(29,15),each=2),Ninds=rep(c(9,7),each=2),Sex=c(0,1),Age=3,Stage=3))

# write InitIndsFile to fil
write.table(init_df_29_15, file=paste0(dirpath,'Inputs/InitInds_29_15.txt'), sep='\t', row.names=F, quote=F)

# Set initialisation
init_29_15 <- Initialise(InitType = 2,       # from loaded species distribution map
                         InitIndsFile = 'InitInds_29_15.txt')



#-----------------
# SIMULATION
#-----------------

# This time, we set a more reasonable number of replicates with 100 runs and simulate range expansion over 100 years. We specify that output should be stored for each year.

# Simulations
RepNb <- 100

sim <- Simulation(Simulation = 0, # ID
                  Replicates = RepNb, # number of replicate runs
                  Years = 100, # number of simulated years
                  OutIntPop = 1, # output interval
                  OutIntOcc = 1,
                  OutIntRange = 1)


#-----------------
# PARAMETER MASTER
#-----------------

# When defining the parameter master objects for our two scenarios, we take care to provide two different batchnum to avoid any overwriting of file output. We set a seed for easy replicability of results such that all of us should obtain the same results.

# RangeShifter parameter master object for single-site reintroduction
s_29 <- RSsim(batchnum = 1, land = land, demog = demo, dispersal = disp, simul = sim, init = init_29, seed = 324135)

# RangeShifter parameter master object for multi-site reintroduction
s_29_15 <- RSsim(batchnum = 3, land = land, demog = demo, dispersal = disp, simul = sim, init = init_29_15, seed = 324135)



#-----------------
# RUN SIMULATION
#-----------------

# Run single-site simulations
RunRS(s_29, dirpath)

# Run multi-site simulations
RunRS(s_29_15, dirpath)


#-----------------
# PLOT RESULTS
#-----------------

#---------------
# Single-site reintroduction
#---------------

# general output of population size + occupancy
par(mfrow=c(1,2))
plotAbundance(s_29,dirpath,sd=T, rep=F)
plotOccupancy(s_29, dirpath, sd=T, rep=F)


# To obtain a deeper understanding of the reintroduction success, we look at different colonisation metrics offered in the function ColonisationStats(). Specifically, we calculate colonisation probability for year 100 (the probability of a patch to be occupied after 100 years) and the time to colonisation.

# Colonisation metrics, the computation of these may take some time
col_stats_29 <- ColonisationStats(s_29, dirpath, years = 100, maps = T)

# We can extract and map occupancy probability.
# mean occupancy probability in year 100
head(col_stats_29$occ_prob)

# map occupancy probability
mycol_occprob <- colorRampPalette(c('blue','orangered','gold'))
levelplot(col_stats_29$map_occ_prob, margin=F, scales=list(draw=FALSE), at=seq(0,1,length=11), col.regions=mycol_occprob(11))

# map colonisation time
mycol_coltime <- colorRampPalette(c('orangered','gold','yellow','PowderBlue','LightSeaGreen'))
levelplot(col_stats_29$map_col_time, margin=F, scales=list(draw=FALSE), at=c(-9,seq(-.001,100,length=11)), col.regions=c('blue',mycol_coltime(11)))


#---------------
# Multi-site reintroduction
#---------------

# general output of population size + occupancy
par(mfrow=c(1,2))
plotAbundance(s_29_15,dirpath,sd=T, rep=F)
plotOccupancy(s_29_15, dirpath, sd=T, rep=F)

# Colonisation metrics
col_stats_29_15 <- ColonisationStats(s_29_15, dirpath, years = 100, maps = T)

# map occupancy probability
levelplot(col_stats_29_15$map_occ_prob, margin=F, scales=list(draw=FALSE), at=seq(0,1,length=11), col.regions=mycol_occprob(11))

# map colonisation time
levelplot(col_stats_29_15$map_col_time, margin=F, scales=list(draw=FALSE), at=c(-9,seq(-.001,100,length=11)), col.regions=c('blue',mycol_coltime(11)))


#-----------------
# EXTINCTION PROBABILITY
#-----------------

# Extinction probability at a specific time can be defined as the proportion of replicate simulation runs without viable population at a specific point in time. We can extract this information from the population output file.
# Accordingly, the mean time to extinction can then be defined as the mean time across all replicates when the population went extinct. Again, this information can be extracted from the population output file.

# For convenience,  we define two own functions (using the construct function()) for calculating extinction probability and mean time to extinction.

# Define a function for calculating extinction probability
Calc_ExtProb <- function(pop_df,s) {
  require(dplyr)
  require(tidyr)
  
  pop_df %>%
    group_by(Rep,Year) %>%
    # Sum individuals over all cells per year and replicate
    summarise(sumPop = sum(NInd), .groups='keep') %>%
    group_by(Year) %>%
    # Average extinction probability (1 minus the proportion of replicates with surviving populations)
    summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb) %>%
    # Make sure that data frame is filled until last year of simulation
    right_join(tibble(Year = seq_len(s@simul@Years)), by='Year') %>% replace_na(list(extProb=1))
}

# Define a function for calculating mean time to extinction
Calc_ExtTime <- function(pop_df) {
  require(dplyr)
  require(tidyr)
  
  pop_df %>%
    group_by(Rep,Year) %>%
    # Sum individuals over all cells per year and replicate    
    summarise(sumPop = sum(NInd), .groups='keep') %>% 
    # Identify in which year they go extinct
    filter(sumPop==0) %>% 
    pull(Year) %>% mean
}

# We can now use these function on our simulation output.

#---------------
# Single-site reintroduction
#---------------

# read population output file into a data frame
pop_29 <- readPop(s_29, dirpath)

# extinction probability
extProb_29 <- Calc_ExtProb(pop_29,s_29)

# Plot extinction probabilities
ggplot(data = extProb_29, mapping = aes(x = Year, y = extProb)) + 
  geom_line() +
  ylim(0,1)

# mean time to extinction
Calc_ExtTime(pop_29)


#---------------
# Multi-site reintroduction
#---------------

# read population output file into a data frame
pop_29_15 <- readPop(s_29_15, dirpath)

# extinction probability
extProb_29_15 <- Calc_ExtProb(pop_29_15,s_29_15)

# mean time to extinction
Calc_ExtTime(pop_29_15)


# Compare extinction probabilities for single- and multi-site reintroduction

# Join extinction probabilities in a single data frame
extProb_scens <- bind_rows(extProb_29 %>% add_column(Scenario = "10 ind. Kintyre "),
                           extProb_29_15 %>% add_column(Scenario = "32 ind. Kintyre + Aberdeenshire"))

ggplot(data = extProb_scens, mapping = aes(x = Year, y = extProb, color=Scenario)) + 
  geom_line(size=2) +
  ylim(0,1)
