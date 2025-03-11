

find_first_linear_chunk <- function(dataframe, gastype, kstar = 0.4, which.chunk = "first", length.min){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[3],]

  mydf <- data.frame(POSIX.time = dataframe$POSIX.time,
                     time = as.numeric(dataframe$POSIX.time-first(dataframe$POSIX.time)),
                     conc = dataframe[[gastype]])

  df.R2 <- NULL
  for(t in seq(length.min, length(mydf$conc))){
    lm.tmp <- lm(formula = conc~time, data = mydf[mydf$time<=t,])
    r2 <- summary(lm.tmp)$adj.r.squared

    df.R2 <- rbind(df.R2, data.frame(t = t,
                                     r2 = r2))
  }

  ind_bests <- which(df.R2$r2>quantile(df.R2$r2, 0.9))
  ind_not_1 <- which(diff(ind_bests)!=1)
  if(length(ind_not_1)==0){
    ind_best <- which.max(df.R2$r2)
  } else {ind_best <- ind_bests[min(ind_not_1)-1]}

  t_selected <- df.R2$t[ind_best]

  plot(mydf$time, mydf$conc)
  points(mydf$time[mydf$time<=t_selected], mydf$conc[mydf$time<=t_selected], col='red')


  return(t_selected)
}
