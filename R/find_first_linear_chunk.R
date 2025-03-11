

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

  ind_bests <- which(df.stats$r2>quantile(df.stats$r2, 0.9))
  ind_not_1 <- which(diff(ind_bests)!=1)
  if(length(ind_not_1)==0){
    ind_best <- which.max(df.stats$r2)
  } else {ind_best <- ind_bests[min(ind_not_1)-1]}

  t_selected <- df.stats$t[ind_best]

  # plot(mydf$time, mydf$conc)
  # points(mydf$time[mydf$time<=t_selected], mydf$conc[mydf$time<=t_selected], col='red')


  first_linear_chunk <- data.frame(UniqueID = unique(dataframe$UniqueID),
                       start = 1,
                       end = t_selected)
  first_linear_chunk$obs.length <- first_linear_chunk$end-first_linear_chunk$start

  first_linear_chunk$start.time <- mydf$POSIX.time[first_linear_chunk$start]


  return(first_linear_chunk)
}




find_first_linear_chunk.loop <- function(x, list_of_dataframes, gastype, length.min) {

  # Function to apply in the loop. Adapt parameters to your needs.
  chunk.out <- find_first_linear_chunk(dataframe = list_of_dataframes[[x]], gastype, length.min)

  return(chunk.out)
}
