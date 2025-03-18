test_that("plot overview works", {

  library(aquaGHG)
  # plot incubations overview
  p <- plot.incubations(dataframe = mydata_all)
  print(p)

  # automatic inspection of CO2 data and flux calculation
  CO2_flux <- automaticflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 0, gastype = "CO2dry_ppm",
                            fluxSeparation = F, displayPlots = T,
                            method = "trust.it.all", force.separation = F) # "trust.it.all" or "focus.on.linear"

})
