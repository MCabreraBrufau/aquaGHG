

# Demo


# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(tidyverse)
library(egg)


# --- install goFlux if needed ---
library(devtools)
# install_github("Qepanna/goFlux")

library(goFlux)


# Sourcing all functions in repo
repo_root <- dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
files.sources = list.files(path = paste0(repo_root,"/R/functions/"), full.names = T)
for (f in files.sources){source(f)}

files.sources = list.files(path = paste0(repo_root,"/R"), full.names = T, pattern = ".R")
for (f in files.sources){source(f)}


# Loading data
setwd(repo_root)
load(file = "data/data_example_3.RData")


# Loading auxfile table
myauxfile = read.csv("data/myauxfile.csv")
myauxfile$start.time <- as.POSIXct(myauxfile$start.time, tz = 'UTC')

mydata_manID <- clickflux(mydata_all = mydata, myauxfile = myauxfile, gastype = "CO2dry_ppm")



# Calculate fluxes
CO2_results_manID <- goFlux(mydata_manID, "CO2dry_ppm")
H2O_results_manID <- goFlux(mydata_manID, "H2O_ppm")
CH4_results_manID <- goFlux(mydata_manID, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res_manID <- best.flux(CO2_results_manID, criteria)
H2O_flux_res_manID <- best.flux(H2O_results_manID, criteria)
CH4_flux_res_manID <- best.flux(CH4_results_manID, criteria)

