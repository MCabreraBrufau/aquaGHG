

clickflux <- function(mydata_all, myauxfile, gastype){

  # ---- Process incubations and compute fluxes ----

  # Define the measurements' window of observation
  # mydata_ow <- obs.win(inputfile = mydata_all, auxfile = myauxfile,
  #                      obs.length = myauxfile$duration, shoulder = 2)

  # Split data into separate dataframes for each incubation to ease following steps
  mydata_ow <- mydata_all %>% group_split(UniqueID) %>% as.list()

  # ----------- Compute fluxes after manual selection of CO2 data

  # Manually identify measurements by clicking on the start and end points
  mydata_manID <- lapply(seq_along(mydata_ow), goFlux::click.peak.loop,
                         flux.unique = mydata_ow,
                         gastype = gastype,
                         plot.lim = c(200,1000)) %>%
    map_df(., ~as.data.frame(.x))


  # Additional auxiliary data required for flux calculation.
  mydata_manID <- mydata_manID %>%
    left_join(myauxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))

  # Add instrument precision for each gas
  mydata_manID <- mydata_manID %>%
    mutate(CO2_prec = first(mydata_all$CO2_prec), CH4_prec = first(mydata_all$CH4_prec),
           N2O_prec = first(mydata_all$N2O_prec), H2O_prec = first(mydata_all$H2O_prec))

  return(mydata_manID)
}





