

find_first_linear_chunk <- function(dataframe, gastype, length.min){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[3],]

  mydf <- data.frame(POSIX.time = dataframe$POSIX.time,
                     time = as.numeric(dataframe$POSIX.time-first(dataframe$POSIX.time)),
                     conc = dataframe[[gastype]])

  df.stats <- NULL
  for(t in seq(length.min, length(mydf$conc))){
    lm.tmp <- lm(formula = conc~time, data = mydf[mydf$time<=t,])
    df.stats <- rbind(df.stats, data.frame(t = t,
                                           r2 = summary(lm.tmp)$adj.r.squared,
                                           pval = summary(lm.tmp)$coefficients[2,4]))
  }

  if(max(df.stats$r2)<0.3){
    message("... no linear chunk could be found, the entire incubation was selected")
    first_linear_chunk <- data.frame(UniqueID = unique(dataframe$UniqueID),
                                     start = 1,
                                     end = length(mydf$conc))
  } else {
    ind_bests <- which(df.stats$r2>quantile(df.stats$r2, 0.9) & df.stats$pval<0.05)
    ind_not_1 <- which(diff(ind_bests)>10)
    if(length(ind_not_1)==0){
      ind_best <- which.max(df.stats$r2)
    } else {ind_best <- ind_bests[min(ind_not_1)-1]}
    t_selected <- df.stats$t[ind_best]

    # p_r2 <- ggplot(df.stats, aes(t, r2))+geom_point()+theme_article()+xlab("time [secs]")+
    #   geom_vline(xintercept = t_selected, colour = "red")+xlim(c(0,t))+ggtitle("LM.r2 over a growing time window")
    # p_conc <- ggplot(mydf, aes(time, conc))+geom_point()+theme_article()+
    #   geom_point(data = mydf[mydf$time<=t_selected,], aes(time, conc), colour = "red")+
    #   ggtitle(unique(dataframe$UniqueID))+geom_vline(xintercept = t_selected, colour = "red")
    # ggarrange(p_conc, p_r2)+xlim(c(0,t))

    first_linear_chunk <- data.frame(UniqueID = unique(dataframe$UniqueID),
                                     start = 1,
                                     end = t_selected)
  }


  first_linear_chunk$obs.length <- first_linear_chunk$end-first_linear_chunk$start
  first_linear_chunk$start.time <- mydf$POSIX.time[first_linear_chunk$start]

  return(first_linear_chunk)
}




find_first_linear_chunk.loop <- function(x, list_of_dataframes, gastype, length.min) {

  # Function to apply in the loop. Adapt parameters to your needs.
  chunk.out <- find_first_linear_chunk(dataframe = list_of_dataframes[[x]], gastype, length.min)

  return(chunk.out)
}
