

# Demo

# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(tidyverse)
library(egg)

repo_root <- dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
files.sources = list.files(path = paste0(repo_root,"/R/functions/"), full.names = T)
for (f in files.sources){source(f)}


setwd(repo_root)
load(file = "data/data_example_3.RData")

flagged <- flux_separator(my_incub = mydata, UniqueID = first(mydata$UniqueID),
                          kstar = 0.6, doPlot = T)


