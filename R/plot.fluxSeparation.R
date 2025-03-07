

plot.fluxSeparation <- function(dataframe, gastype, kstar = 0.5) {

  ## Check dataframe ####
  if(missing(dataframe)) stop("'dataframe' is required")
  if(!is.null(dataframe) & !is.data.frame(dataframe)){
    stop("'dataframe' must be of class data.frame")}

  ### gastype and match in dataframe ####
  if(missing(gastype)) stop("'gastype' is required")
  if(!is.null(gastype) & !is.character(gastype)) stop("'gastype' must be a character string")
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""),
                c("CO2dry_ppm", "COdry_ppb", "CH4dry_ppb", "N2Odry_ppb", "NH3dry_ppb", "H2O_ppm")))){
    stop("'gastype' must be of class character and one of the following: 'CO2dry_ppm', 'COdry_ppm', 'CH4dry_ppb', 'N2Odry_ppb', 'NH3dry_ppb' or 'H2O_ppm'")}
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe)))){
    stop("'dataframe' must contain a column that matches 'gastype'")}
  if(any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe))) &
     !is.numeric(dataframe[,gastype][[1]])){
    stop("The column that matches 'gastype' in 'dataframe' must be of class numeric")}

  ## Check kstar ####
  if(is.null(kstar)) stop("'kstar' is required") else{
    if(!is.numeric(kstar)) stop("'kstar' must be of class numeric") else{
      if(kstar < 0) stop("'kstar' cannot be a negative value")
      if(kstar > 1) stop("'kstar' cannot be > 1")}}

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

  # seq.rep: A wrapper function that merges seq and rep, to repeat a sequence
  # or sequence a repeat.
  seq.rep <- function(from, by, n.rep, length.seq, rep.seq = F) {

    if (rep.seq == F) {

      sequence <- seq(from = from, by = by, length.out = length.seq)
      out.ls <- list()
      for (i in 1:length(sequence)) {
        out.ls[[i]] <- rep(sequence[i], n.rep) }
      unlist(out.ls)

    } else {

      repetition <- rep(from, n.rep)
      out.ls <- list()
      for (i in 1:length(repetition)) {
        out.ls[[i]] <- seq(from = repetition[i], by = by, length.out = length.seq) }
      unlist(out.ls)

    }
  }


  # FUNCTION STARTS ####


  # Define gas units
  if(gastype == "CO2dry_ppm") gas.unit <- "ppm"
  if(gastype == "CH4dry_ppb") gas.unit <- "ppb"
  if(gastype == "N2Odry_ppb") gas.unit <- "ppb"
  if(gastype == "COdry_ppb") gas.unit <- "ppb"
  if(gastype == "NH3dry_ppb") gas.unit <- "ppb"
  if(gastype == "H2O_ppm") gas.unit <- "ppm"

  # Define y axis legend on plots
  if(gastype == "CO2dry_ppm") ylab <- ylab(expression(bold(CO["2"]*" dry (ppm)")))
  if(gastype == "CH4dry_ppb") ylab <- ylab(expression(CH["4"]*" dry (ppb)"))
  if(gastype == "N2Odry_ppb") ylab <- ylab(expression(N["2"]*"O dry (ppb)"))
  if(gastype == "COdry_ppb") ylab <- ylab(expression(CO*" dry (ppb)"))
  if(gastype == "NH3dry_ppb") ylab <- ylab(expression(NH["3"]*" dry (ppb)"))
  if(gastype == "H2O_ppm") ylab <- ylab(expression(H["2"]*"O (ppm)"))


  # Initializing column with flag for diffusion
  dataframe$flag_diffusion <- TRUE

  # Create a list of dataframe (by UniqueID)
  data_split <- dataframe %>%
    group_by(UniqueID) %>%
    group_split()


  # Loop through list of data frames (by UniqueID)
  pboptions(char = "=")
  plot_list <- pblapply(seq_along(data_split), function(f) {

    # flag diffusion
    data_split[[f]]$flag_diffusion <- flag_diffusion(dataframe = data_split[[f]], kstar = kstar)

    # Plot limits
    xmax <- max(na.omit(data_split[[f]]$Etime))
    xmin <- min(na.omit(data_split[[f]]$Etime))
    xdiff <- xmax - xmin

    ymax <- max(na.omit(data_split[[f]][, gastype]))
    ymin <- min(na.omit(data_split[[f]][, gastype]))
    ydiff <- ymax - ymin

    # Content of plot
    Etime <- data_split[[f]]$Etime
    gas_meas <- Reduce("c", data_split[[f]][, gastype])
    flag <- data_split[[f]]$flag_diffusion
    component <- ""
    component[flag] <- "diffusion"
    component[!flag] <- "ebullition"
    plot_data <- cbind.data.frame(gas_meas, Etime, flag, component)

    incubationID <- unique(data_split[[f]]$UniqueID)

    # Draw plot ####
    plot <- ggplot(plot_data, aes(x = Etime)) +
      geom_point(aes(y = gas_meas, colour = component)) +
      scale_color_manual(values = c("grey20", "#FF7F50")) +

      # Make the plot pretty
      xlab("Time (sec)") + ylab +
      scale_x_continuous(breaks = seq(-60, max(Etime), 60),
                         minor_breaks = seq(-60, max(Etime)+60, 10)) +
      theme_bw() +
      theme(axis.title.x = element_text(size = 10, face = "bold"),
            axis.title.y = element_text(size = 10, face = "bold"))+
      ggtitle(incubationID)

    return(plot)
  })
}
