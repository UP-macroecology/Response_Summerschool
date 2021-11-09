# Example SDM workflow
# RTG Response summer school 2021
# (c) Damaris Zurell, Univ. Potsdam
# More elaborate materials: https://damariszurell.github.io/SDM-Intro/ and https://damariszurell.github.io/EEC-MGC/ 

#--------------------------------------------------------------------------------
#
#         Load packages
# 
#--------------------------------------------------------------------------------

library(raster)   # Manipulate geographic data
library(corrplot)   # Graphical displays of correlation matrices

# Installing package "mecofun" with helper functions (meant just for teaching purposes)
library(devtools)
devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")
library(mecofun)


#--------------------------------------------------------------------------------
#
#         SET WORKING DRECTORY
# 
#--------------------------------------------------------------------------------

# Set your working directory to the workshop folder
setwd('data/Lehre/workshops_tutorials/Response_2021/Response_Summerschool/SDM')



#--------------------------------------------------------------------------------
#
#         TYPICAL SDM DATA
# 
#--------------------------------------------------------------------------------

# DATA

# Data stem from Citizen Science Atlas UK (https://doi.org/10.1111/geb.12906). 
# Data set contains presence-absence data of the Ring Ouzel in UK for breeding period 2008-2011. Climate data stem from worldclim. 

# Read in presence-absence data
sp_dat <- read.table('data/ATLAS_RingOuzel.txt',header=T)

# Inspect data
summary(sp_dat)

# Map data within Britsh Isle
# Read in background data defining British Isle land mass:
bg <- raster('data/UK_mask.grd')  
# Plot GB land mass:
plot(bg,col='grey',axes=F,legend=F)   
# Plot presences in red and absences in black:
plot(extend(rasterFromXYZ(sp_dat[,1:3]),bg), col=c('black','red'), legend=F,add=T)


#--------------------------------------------------------------------------------
#
#         SIMPLE SDM FITTING (GLM)
# 
#--------------------------------------------------------------------------------


#---------------
# GLM FORMULAS
#---------------

# We first fit a GLM for the bio11 variable assuming a linear relationship:
m1 <- glm(Turdus_torquatus ~ bio11, family="binomial", data= sp_dat)

# We can get a summary of the model:
summary(m1) 

# Some options for more complex model specifications
# Fit a quadratic relationship with bio11:
m1_q <- glm(Turdus_torquatus ~ bio11 + I(bio11^2), family="binomial", data= sp_dat)
summary(m1_q)

# Or use the poly() function:
summary( glm(Turdus_torquatus ~ poly(bio11,2) , family="binomial", data= sp_dat) )

# Fit two variables with second-order polynomials:
summary( glm(Turdus_torquatus ~ poly(bio11,2) + poly(bio8,2), family="binomial", data= sp_dat) )


#---------------
# COLLINEARITY
#---------------

# We first estimate a correlation matrix from the predictors. 
# We use Spearman rank correlation coefficient, as we do not know 
# whether all variables are normally distributed.
cor_mat <- cor(sp_dat[,-c(1:3)], method='spearman')

# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Use select07 method to dentify all pairs of variables that have correlation |rho|>0.7 and remove the less important variable
# Function described in Dormann et al. (2013): http://dx.doi.org/10.1111/j.1600-0587.2012.07348.x

# Run select07()
var_sel <- select07(X=sp_dat[,-c(1:3)], 
                    y=sp_dat$Turdus_torquatus, 
                    threshold=0.7)

# Check out the structure of the resulting object:
str(var_sel)

# We extract the names of the weakly correlated predictors ordered by the univariate variable importance in terms of AIC:
pred_sel <- var_sel$pred_sel

# How many presence points do we have? Rule of thumb: you need 10 presences per parameter in the model
sum(sp_dat$Turdus_torquatus)


#----------------
# MODEL SELECTION
#----------------

# Fit the full model with 4 parameters:
m_full <- glm( Turdus_torquatus ~ bio11 + I(bio11^2) + bio8 + I(bio8^2), 
               family='binomial', data=sp_dat)

# Inspect the model:
summary(m_full)

# Explained deviance:
expl_deviance(obs = sp_dat$Turdus_torquatus,
              pred = m_full$fitted)
