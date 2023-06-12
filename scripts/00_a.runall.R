# Description  ------------------------------------------------------------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Run all script 
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created January 27, 2023
# Why: 
  # Run all scripts in order to 
    # produce outputs, 
    # intermediate data, 
    # figures, and tables. 
# Requires 
  # Necessary scripts 
# NOTES: 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Remove all stored objects and garbage collector  ------------------------
rm(list = ls()); gc(verbose = TRUE, reset = FALSE, full = TRUE)

# Set site to check -------------------------------------------------------
# Select the site for which the data will be computed 
site.check = "El Garrapatero" 
ext.file = "_EG"


# Run each script in order ------------------------------------------------
source("scripts/0.0_initialize.R")                       # Libraries and functions

source("scripts/fitness_surfaces/fit.land.analysis.R")   # Further preparation 
source("scripts/fitness_surfaces/fit.land.analysis_2.R") # This step is long
source("scripts/fitness_surfaces/fit.land.analysis_3.R") # Get fitness landscape and prospective selection 

source("scripts/fitness_surfaces/plot_fitness_no_model.R") # Get fitness landscape without model 
source("scripts/adaptive_landscapes/adaptive.landscape.calculation.R") # Get adaptive landscape COMBINED with fitness landscape and with no model. Will be long if you estimate the adaptive landscape again! 900 iterations will be done
source("scripts/fitness_surfaces/fit.land.walk.R")       # Walking the fitness landscape to get information regarding peak distance and valleys 

source("scripts/climate/climate.analysis.R")             # Climate analysis 
source("scripts/maps/Maps sites Galapagos drawing net and sites.R")    # Maps for sites and islands
source("scripts/CMR/cmr.test.R")                         # Capture-mark-recapture analysis (takes a long time)
source("scripts/Repeat.R")                               # Repeatability (takes a long time)
