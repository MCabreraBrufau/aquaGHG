

plot.incubations <- function(dataframe) {

  ## Check dataframe ####
  if(missing(dataframe)) stop("'dataframe' is required")
  if(!is.null(dataframe) & !is.data.frame(dataframe)){
    stop("'dataframe' must be of class data.frame")}

  ### Etime  ####
  if(!any(grepl("\\<Etime\\>", names(dataframe)))) stop("'dataframe' must contain 'Etime'")
  if(any(grepl("\\<Etime\\>", names(dataframe))) & !is.numeric(dataframe$Etime)){
    stop("'Etime' in 'dataframe' must be of class numeric (or integer)")}


  ### POSIX.time ####
  if(!any(grepl("\\<POSIX.time\\>", names(dataframe)))) stop("'dataframe' must contain 'POSIX.time'")
  if(any(grepl("\\<POSIX.time\\>", names(dataframe))) & !is.POSIXct(dataframe$POSIX.time)){
    stop("'POSIX.time' in 'dataframe' must be of class POSIXct")}

  ### UniqueID (or chamID) ####
  if(!any(grepl(paste(c("\\<UniqueID\\>", "\\<chamID\\>"), collapse = "|"), names(dataframe)))){
    stop("'dataframe' must contain 'UniqueID'")}

  # FUNCTION STARTS ####

  # Define gas units
  # if(gastype == "CO2dry_ppm") gas.unit <- "ppm"
  # if(gastype == "CH4dry_ppb") gas.unit <- "ppb"
  # if(gastype == "N2Odry_ppb") gas.unit <- "ppb"
  # if(gastype == "COdry_ppb") gas.unit <- "ppb"
  # if(gastype == "NH3dry_ppb") gas.unit <- "ppb"
  # if(gastype == "H2O_ppm") gas.unit <- "ppm"

  # find out which variables are present in dataframe
  vars <- c("CO2dry_ppm", "CH4dry_ppb", "N2Odry_ppb", "COdry_ppb", "NH3dry_ppb", "H2O_ppm")
  vars_in <- NULL
  for (var in vars){
    if(any(grepl(paste0("\\<",var,"\\>"), names(dataframe)))){vars_in <- c(vars_in, var)}
  }

  dataframe_gath <- gather(dataframe[,c("UniqueID","Etime",vars_in)], variable, value, -Etime, -UniqueID)

  # Create a list of dataframe (by UniqueID)
  data_split <- dataframe_gath %>%
    group_by(UniqueID) %>%
    group_split()


  # Loop through list of data frames (by UniqueID)
  pboptions(char = "=")
  plot_list <- pblapply(seq_along(data_split), function(f) {

    # Content of plot
    Etime <- data_split[[f]]$Etime
    gas_meas <- data_split[[f]]$value
    variable <- data_split[[f]]$variable

    plot_data <- cbind.data.frame(gas_meas, Etime, variable)

    incubationID <- unique(data_split[[f]]$UniqueID)

    # Draw plot ####
    plot <- ggplot(plot_data, aes(x = Etime, y = gas_meas)) +
      geom_point(size=0.5)+
      geom_path()+
      facet_wrap(.~variable, scales = "free_y", ncol = 1)+

      # Make the plot pretty
      xlab("Time (sec)") +
      ylab("Measured gas") +
      scale_x_continuous(breaks = seq(-60, max(Etime), 60),
                         minor_breaks = seq(-60, max(Etime)+60, 10)) +
      theme_article() +
      theme(axis.title.x = element_text(size = 10, face = "bold"),
            axis.title.y = element_text(size = 10, face = "bold"))+
      ggtitle(incubationID)

    return(plot)
  })
}
