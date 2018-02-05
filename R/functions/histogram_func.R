# functions to create histograms for spont and other conditions
# generate 500 rows x 100 of histogram values (x1...x100) 
# of values mean normalised by baseline
spontaneousHist <- function(spont_times=NULL, sweeps=100, iteration=500, binwidth=0.01,
                            pre = -0.5, post = 1){
  max_time <- max(spont_times) - post
  
  nbins <- length(seq(pre,post,binwidth)) - 1
  hist_df <- as.data.frame(matrix(nrow = iteration, ncol = nbins))
  for(i in 1:iteration){
    spike_times <- numeric(0)
    for(j in 1:sweeps){
      stim_time <- sample(seq(0, max_time, 0.001), size = 1)
      spikes_in_sweep <- spont_times[spont_times >= stim_time + pre & spont_times <= stim_time + post]
      spikes_in_sweep <- spikes_in_sweep - stim_time
      spike_times <- c(spike_times, spikes_in_sweep)
    }
    histo <- hist(spike_times, breaks = seq(pre, post, binwidth), plot = F)
    baseline_avg <- mean(histo$counts[1:length(histo$mids[histo$mids < 0])])
    counts <- (histo$counts / baseline_avg) * 100
    hist_df[i, ] <- counts
  }
  
  dfs <- list(full=hist_df, post=hist_df[, (length(histo$mids[histo$mids < 0])+1):(nbins)], mids = histo$mids)
  
  return(dfs)
}

# function for each stimulus condition, histogram values normalised to 100% of baseline
stimulusHistogram <- function(spike_times = NULL, binwidth=0.01, pre=-0.5, post=1){
  spike_times <- spike_times[spike_times > pre & spike_times < post]
  nbins <- length(seq(pre, post, binwidth)) - 1
  histo <- hist(spike_times, breaks=seq(pre, post, binwidth), plot = F)
  baseline_avg <- mean(histo$counts[1:length(histo$mids[histo$mids < 0])])
  
  counts <- (histo$counts / baseline_avg) * 100
  
  dfs <- list(full=counts, post=counts[(length(histo$mids[histo$mids < 0])+1):(nbins)])
  
  return(dfs)
}

