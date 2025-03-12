


find_linear_chunk <- function(dataframe, gastype, kstar = 0.4, which.chunk = "first", length.min){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[2],]

  # computing density probability of first derivative
  d_df <- get_dCdt_density(dataframe, gastype)
  d <- d_df[[1]]
  mydf <- d_df[[2]]

  half_dens_max <- kstar*max(d$y)
  ind_over_half_dens_max <- which(d$y>half_dens_max)

  # we define lower and upper boundaries of the "main slope" as the slopes with
  # a density of probability higher than half of the maximum density
  lower_bound <- d$x[first(ind_over_half_dens_max)]
  upper_bound <- d$x[last(ind_over_half_dens_max)]

  mydf$is_linear <- F
  mydf$is_linear[mydf$dydt>lower_bound & mydf$dydt<upper_bound] <- T

  vect <- mydf$is_linear
  jumps <- NULL
  if(vect[1]==T){jumps <- c(1)}
  for(i in seq(2,length(vect))){
    if(vect[i] != vect[i-1]){
      jumps <- c(jumps, i)
    }
  }
  if(last(vect)==T){jumps <- c(jumps, i)}


  chunks <- data.frame(UniqueID = unique(dataframe$UniqueID),
                       start = jumps[seq(1,length(jumps),2)],
                       end = jumps[seq(2,length(jumps),2)])
  chunks$obs.length <- chunks$end-chunks$start

  chunks$start.time <- mydf$POSIX.time[chunks$start]

  chunks <- chunks[chunks$obs.length>=length.min,]

  if (which.chunk == "first"){chunk.out <- chunks[1,]}
  if (which.chunk == "last"){chunk.out <- chunks[length(chunks$start),]}
  if (which.chunk == "longest"){chunk.out <- chunks[which.max(chunks$obs.length),]}

  return(chunk.out)
}


find_linear_chunk.loop <- function(x, list_of_dataframes, gastype, kstar = 0.4, which.chunk, length.min) {

  # Function to apply in the loop. Adapt parameters to your needs.
  chunk.out <- find_linear_chunk(dataframe = list_of_dataframes[[x]], gastype, kstar, which.chunk, length.min)

  return(chunk.out)
}
