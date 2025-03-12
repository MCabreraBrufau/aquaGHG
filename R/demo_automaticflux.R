

# Demo


# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(tidyverse)
library(egg)
library(goFlux)
library(devtools)
library(zoo)
library(pbapply)
library(ggnewscale)

source(file = "R/automaticflux.R")
source(file = "R/get_dxdy.R")
source(file = "R/join_auxfile_with_data.loop.R")
source(file = "R/plot.incubations.R")
source(file = "R/find_first_linear_chunk.R")
source(file = "R/plot.fluxSeparation.R")
source(file = "R/flag_diffusion.R")
source(file = "R/get_dCdt_density.R")
source(file = "R/flux.separator.R")
source(file = "R/separated.flux.plot.R")

# Loading data
mydata_all <- NULL
fs <- list.files(path = "data/",pattern = ".RData", full.names = T)
for(f in fs){
  load(file = f)
  mydata$Etime <- as.numeric(mydata$Etime)
  mydata_all <- rbind(mydata_all, mydata)
  rm(mydata)
}

# Loading auxfile table
myauxfile = read.csv("data/myauxfile.csv")
myauxfile$start.time <- as.POSIXct(myauxfile$start.time, tz = 'UTC', format="%d/%m/%Y %H:%M")


# plot incubations overview
p <- plot.incubations(mydata_all)
print(p)
# to save these plots in a dedicated path, do
# gg_save_pdf(list = p, path = , filename = "myfilename.pdf")


# automatic inspection of CO2 data and flux calculation
CO2_flux <- automaticflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 30, gastype = "CO2dry_ppm",
                          fluxSeparation = F, displayPlots = T,
                          method = "focus.on.linear") # trust.it.all or focus.on.linear



# take a look at CH4 data and identify possible ebullition events
p <- plot.fluxSeparation(dataframe = mydata_all, gastype = "CH4dry_ppb", kstar = 0.4)
print(p)


# automatic inspection of CO2 data and flux calculation
CH4_flux <- automaticflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 30, gastype = "CH4dry_ppb",
                          fluxSeparation = T, displayPlots = T,
                          method = "trust.it.all") # trust.it.all or focus.on.linear




