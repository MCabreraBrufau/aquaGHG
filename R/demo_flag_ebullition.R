

# Demo

# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(tidyverse)
library(egg)

source("R/flag_diffusion.R")
source("R/plot.fluxSeparation.R")
source("R/get_dxdy.R")

load(file = "data/data_example_2.RData")
mydata$Etime <- as.numeric(mydata$Etime)

# take a look at CH4 data and identify possible ebullition events
# user can play with kstar to be more or less permissive.
# Note that kstar must be chosen between 0.1 and 0.9
# kstar close to 0 is very permissive, close to 1 is very strict.
p <- plot.fluxSeparation(dataframe = mydata, gastype = "CH4dry_ppb", kstar = 0.2)
print(p)

# get flag column in the dataframe mydata
mydata$flag_diffusion <- flag_diffusion(dataframe = mydata, kstar = 0.1)


