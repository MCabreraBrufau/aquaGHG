test_that("plot overview works", {

  library(aquaGHG)
  library(lubridate)
  library(tidyr)
  library(dplyr)
  library(pbapply)
  library(ggplot2)
  library(egg)

  # plot incubations overview
  p <- plot_incubations(dataframe = mydata_all)
  print(p)
})



test_that("automatic CO2 calculation works with method 'trust.it.all' ", {

  library(aquaGHG)
  library(lubridate)
  library(tidyr)
  library(dplyr)
  library(pbapply)
  library(ggplot2)
  library(egg)
  library(goFlux)
  library(purrr)


  # automatic inspection of CO2 data and flux calculation
  CO2_flux <- automaticflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 30, gastype = "CO2dry_ppm",
                            fluxSeparation = FALSE, displayPlots = TRUE,
                            method = "trust.it.all", force.separation = FALSE) # "trust.it.all" or "focus.on.linear"

})




test_that("automatic CO2 calculation works with method 'trust.it.all' ", {

  library(aquaGHG)
  library(lubridate)
  library(tidyr)
  library(dplyr)
  library(pbapply)
  library(ggplot2)
  library(egg)
  library(goFlux)
  library(purrr)


  # automatic inspection of CO2 data and flux calculation
  CO2_flux <- automaticflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 30, gastype = "CO2dry_ppm",
                            fluxSeparation = FALSE, displayPlots = TRUE,
                            method = "focus.on.linear", force.separation = FALSE) # "trust.it.all" or "focus.on.linear"

})





test_that("automatic CH4 calculation works with method 'trust.it.all' and flux separation ", {

  library(aquaGHG)
  library(lubridate)
  library(tidyr)
  library(dplyr)
  library(pbapply)
  library(ggplot2)
  library(egg)
  library(goFlux)
  library(purrr)
  library(ggnewscale)


  # automatic inspection of CO2 data and flux calculation
  CH4_flux <- automaticflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 30, gastype = "CH4dry_ppb",
                            fluxSeparation = T, displayPlots = TRUE,
                            method = "trust.it.all", force.separation = FALSE) # "trust.it.all" or "focus.on.linear"

})


