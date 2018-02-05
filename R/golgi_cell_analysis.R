# Analysis of golgi cells:
#   1. - Spontaneous: Hz, ISIs, CV of log ISI, can be classified as golgi
#   2. - responses signif to stim. count of cells that are to both
#   3. - difference in timing, signif difference between onset/offset times by stim
#   4. - depression to MCP is longer than ISI, individual boxplot - kruskal wallis signif
#   5. - is there a signif difference in those where the MCP stim has no signif effect (control)
#   6. - short latency activations, timing in response to signif 1MCP, averages, max, min etc.


# libraries needed
library(dplyr)
library(ggplot2)
library(reader)
library(cpm)


 
# standard error function 
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# functions to deterimine isi from stim spikes
stimISI <- function(stim_times = NULL, n_pulses = NULL, driven=NULL){
  
  if(driven){
    # any spikes during stim period
    if(any(stim_times > 0 & stim_times < n_pulses*0.0033)){
      final_spike <- last(stim_times[stim_times > 0 & stim_times < n_pulses*0.0033])
      spikes_after <- stim_times[stim_times > n_pulses*0.0033]
      # if there are no after that then can't return an ISI
      if(length(spikes_after)==0){
        return(NULL)
      } else {
        # return the difference between the first spike and the next
        return(first(spikes_after) - final_spike)
      }
    } else {
      return(NULL)
    }
  } else {
    if(any(stim_times > 0 & stim_times < n_pulses*0.0033)){
      return(NULL)
    } else {
      final_spike <- last(stim_times[stim_times < 0])
      spikes_after <- stim_times[stim_times > 0]
      if(is.na(final_spike) | length(spikes_after) == 0){
        return(NULL)
      } else {
        return(first(spikes_after) - final_spike)
      }
    }
  }
}

# variables for histogram
pre <- -0.5
post <- 2
binwidth <- 0.01

# set working directory
setwd(choose.dir())

# list all folders
folders <- list.dirs()
folders <- gsub("./", x = folders, replacement = "", fixed = T)
folders <- folders[folders!="."]
cell_num <- length(folders)

# spontanoues data frame
spont_analysis <- data.frame(cell= folders, Hz=rep(NA, cell_num), avg_ISI = rep(NA, cell_num), 
                             CVlogISI = rep(NA, cell_num), depression_to_stim = rep(NA, cell_num))

# significant response to stim
response_analysis <- data.frame(cell = folders, hindlimb = rep(NA, cell_num), forelimb = rep(NA, cell_num),
                                vib = rep(NA, cell_num), MCP1 = rep(NA, cell_num), MCP5 = rep(NA, cell_num))
# onset data frame
onset_df <- data.frame(cell = folders, hindlimb = rep(NA, cell_num), forelimb = rep(NA, cell_num),
                       vib = rep(NA, cell_num), MCP1 = rep(NA, cell_num), MCP5 = rep(NA, cell_num))

# offset dataframe
offset_df <- data.frame(cell = folders, hindlimb = rep(NA, cell_num), forelimb = rep(NA, cell_num),
                        vib = rep(NA, cell_num), MCP1 = rep(NA, cell_num), MCP5 = rep(NA, cell_num))

cell_isis <- as.data.frame(matrix(nrow=0, ncol=5))
colnames(cell_isis) <- c("cell", "condition", "isis", "norm.isi", "signif")

# counts of numbers that have significantly different isis for those that are significant or not
kruskal_isi <- data.frame(cell = folders, response_to_mcp=rep(NA, cell_num), p_value=rep(NA, cell_num))

# spikes in to response 1 MCP
driven_spikes <- as.data.frame(matrix(nrow=0, ncol=2))
colnames(driven_spikes) <- c("cell", "times")

# probability of spikes in 5MCP
driven_probability <- as.data.frame(matrix(nrow=0, ncol=3))
colnames(driven_probability) <- c("cell", "stim_number", "probability")

# cell changepoints
cell_changepoints <- as.data.frame(matrix(nrow=0, ncol=3))
colnames(cell_changepoints) <- c("cell", "condition", "changepoints")

# iterate through each cell folder
for(neuron in folders){
  # set the wd as that cell
  setwd(neuron)
  cat("working on", which(folders==neuron), "out of", length(folders), "\n")
  
  # list the files in the neuron folder
  files <- list.files()
  
  dir.create("general plots")
  
  # spontaneous spike times
  total_file_time <- as.numeric(readLines(grep("spontaneous", files, value = T))[5])
  spont_times <- as.numeric(readLines(grep("spontaneous", files, value = T))[6:file.nrow(grep("spontaneous", files, value = T))])
  
  # stats for spont; Hz
  spont_analysis$Hz[spont_analysis$cell == neuron] <- 1/(total_file_time/length(spont_times))
  
  # stats for spont; ISIs, CV log ISI
  shifted_spont <- c(spont_times[2:length(spont_times)], NA)
  ISIs <- na.omit(shifted_spont - spont_times)
  spont_analysis$avg_ISI[spont_analysis$cell == neuron] <- mean(ISIs)
  # LCV = coefficient of variation for log ISI (in ms) 
  logISI <- log(ISIs*1000)
  LCV <- sd(logISI)/mean(logISI)
  spont_analysis$CVlogISI[spont_analysis$cell == neuron] <- LCV
  
  # plot of ISIs
  spontISI <- ggplot() + geom_histogram(aes(x=ISIs), binwidth = 0.01) + 
    theme_classic() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(limits = c(0,0.5), 
                                                                             labels=function(x)x*1000,
                                                                             name = "interspike interval (ms)",
                                                                             expand = c(0,0))
  dir.create("spont isi plot")
  setwd("spont isi plot")
  ggsave(spontISI, filename = paste0(neuron, " spontaneous ISI plot.pdf"), width = 6, height = 4)
  setwd("..")
  
  # bind the isis to cell isis df
  cell_isis <- rbind(cell_isis, data.frame(cell=rep(neuron, length(ISIs)), condition=rep("spontaneous", length(ISIs)), isis = ISIs,
                                           norm.isi = ISIs/mean(ISIs), signif=rep(NA, length(ISIs))))
  
  # 500 histograms (100 sweeps) of spontaneous
  spontaneous_histo <- spontaneousHist(spont_times = spont_times, pre = pre, post = post)
  smooth_spont <- as.data.frame(matrix(nrow=length(spontaneous_histo$mids)*500, ncol=6))
  colnames(smooth_spont) <- c("time", "value", "value2", "group", "alpha", "condition")
  smooth_spont$time <- spontaneous_histo$mids
  smooth_spont$alpha <- 0.1
  smooth_spont$condition <- "spont"
  
  # iterate through and create a smooth spline for each spontaneous
  start_point <- 1
  end_point <- length(spontaneous_histo$mids)
  length_spont <- end_point
  for(i in 1:nrow(spontaneous_histo$full)){
    smooth_spont$value[start_point:end_point] <- smooth.spline(spontaneous_histo$mids, y=spontaneous_histo$full[i,], spar=0.5)$y
    smooth_spont$value2[start_point:end_point] <- smooth.spline(spontaneous_histo$mids, y=spontaneous_histo$full[i,], spar=0.3)$y
    smooth_spont$group[start_point:end_point] <- i
    start_point <- start_point + length_spont
    end_point <- end_point + length_spont
  }
  
  # generate empty stimulus
  stimulus_histograms <- list()
  full_histograms <- list()

  # stimulus conditions
  colour_key <- list(forelimb="#ff0000", hindlimb="#1d0eba", vib="#23c920", '1MCP'="#8a10de", '5MCP'="#f4a300")
  conditions <- c("hindlimb", "forelimb", "vib", "1MCP", "5MCP")
  for(cond in conditions){
    
    # generate histograms for stimulus conditions if they exist
    if(any(grepl(cond, files))){
      # read in the spike times
      spike_table <- read.table(grep(cond, files, value = T), header = T, sep = "\t",
                                   colClasses = c("numeric", "numeric", "character"), fill = T)
      spike_times <- as.numeric(unlist(strsplit(spike_table[, 3], split = ",")))
      
      # if MCP, split out driven spike and not driven ISIs
      if(grepl("MCP", cond)){
        spike_isi <- strsplit(spike_table[, 3], split = ",")
        spike_isi <- lapply(spike_isi, as.numeric)
        
        pulses <- ifelse(cond=="5MCP", 5, 1)
        driven <- unlist(lapply(spike_isi, stimISI, n_pulses = pulses, driven = T))
        not_driven <- unlist(lapply(spike_isi, stimISI, n_pulses = pulses, driven = F))
        
        if(length(driven)!=0){
          cell_isis <- rbind(cell_isis, data.frame(cell=rep(neuron, length(driven)), condition=rep(paste0(cond, "-excitation"), length(driven)), isis = driven,
                                                   norm.isi = driven/mean(ISIs), signif=rep(NA, length(driven))))
        }
        if(length(not_driven)!=0){
          cell_isis <- rbind(cell_isis, data.frame(cell=rep(neuron, length(not_driven)), condition=rep(paste0(cond, "-no excitation"), length(not_driven)), isis = not_driven,
                                                   norm.isi = not_driven/mean(ISIs), signif=rep(NA, length(not_driven))))
        }
      }
      
      # process the spike times
      stim_hist <- stimulusHistogram(spike_times = spike_times, pre = pre, post = post)
      
      # place it in the stimulus list
      stimulus_histograms[[cond]] <- stim_hist$post
      full_histograms[[cond]] <- stim_hist$full
      
      # smoothing parameter, for MCP stim it smooths too much over time 0
      spar <- ifelse(grepl("MCP", cond), 0.3, 0.5)
      if(grepl("MCP", cond)){
        value2 <- smooth.spline(x=spontaneous_histo$mids, y=stim_hist$full, spar=spar)$y
        value <- rep(NA, length_spont)
      } else {
        value <- smooth.spline(x=spontaneous_histo$mids, y=stim_hist$full, spar=spar)$y
        value2 <- rep(NA, length_spont)
      }
      
      # smooth hind limb
      smooth_stim <- data.frame(time=spontaneous_histo$mids, 
                                value=value,
                                value2=value2,
                                group=rep(length_spont+1, length_spont),
                                alpha=rep(0.8, length_spont), condition=rep(cond, length_spont))
      
      # in case any values go below 0
      smooth_stim$value[smooth_stim$value < 0] <- 0
      smooth_stim$value2[smooth_stim$value2 < 0] <- 0
      
      # bind it spontaneous
      smooth_stim <- rbind(filter(smooth_spont, group%in%sample(1:length_spont, 25)), smooth_stim)
      
      # max for y
      y_max <- ifelse(max(smooth_stim$value) < 200, 200, NA)
      
      
      if(grepl("MCP", cond)){
        stim_plot <- ggplot(smooth_stim, aes(x=time, y=value2, group=group, colour=condition)) +
          geom_line(aes(alpha=alpha), size=1) + scale_alpha_identity() + theme_classic() +
          scale_x_continuous(name="time relative to stimulus (ms)", limits=c(-0.3, 1.8), labels=function(x)x*1000,
                             expand = c(0,0)) + 
          scale_y_continuous(name="% activity compared to baseline", limits=c(0,y_max), expand=c(0,0)) + 
          scale_color_manual(values = c(colour_key[[cond]], "#000000"), guide=guide_legend(title = ""), labels=c(cond, "spontaneous"))
      } else if(cond=="vib"){
        stim_plot <- ggplot(smooth_stim, aes(x=time, y=value, group=group, colour=condition)) +
          geom_line(aes(alpha=alpha), size=1) + scale_alpha_identity() + theme_classic() +
          scale_x_continuous(name="time relative to stimulus (ms)", limits=c(-0.3, 1.8), labels=function(x)x*1000,
                             expand = c(0,0)) + 
          scale_y_continuous(name="% activity compared to baseline", limits=c(0,y_max), expand=c(0,0)) + 
          scale_color_manual(values = c("#000000", colour_key[[cond]]), guide=guide_legend(title = ""), labels=c("spontaneous", cond))
      
      } else {
        stim_plot <- ggplot(smooth_stim, aes(x=time, y=value, group=group, colour=condition)) +
          geom_line(aes(alpha=alpha), size=1) + scale_alpha_identity() + theme_classic() +
          scale_x_continuous(name="time relative to stimulus (ms)", limits=c(-0.3, 1.8), labels=function(x)x*1000,
                             expand = c(0,0)) + 
          scale_y_continuous(name="% activity compared to baseline", limits=c(0,y_max), expand=c(0,0)) + 
          scale_color_manual(values = c(colour_key[[cond]], "#000000"), guide=guide_legend(title = ""), labels=c(cond, "spontaneous"))
      }
      
      # set the wd in the general plot folder
      setwd("general plots")
      ggsave(stim_plot, filename = paste0(cond, " plot.pdf"))
      setwd("..")                      
    }
  }
  
  
  # significance of histogram against spontaneous multivariate normal distribution
  # does the X1...Xn for the stim lie outside a confidence interval given as 1-alpha
  params <- normalParams(spontaneous_histo$post)
  signif_table <- gaussianConfInt(x = stimulus_histograms, mu = params$mu,
                                  Sigma = params$sigma, alpha = 0.0001)
  
  
  
  
  # add the significance to the table
  for(conds in signif_table$condition){
    if(conds=="1MCP"){
      condition <- "MCP1"
    } else if(conds=="5MCP"){
      condition <- "MCP5"
    } else {
      condition <- conds
    }
    response_analysis[response_analysis$cell==neuron, condition] <- signif_table$signif[signif_table$condition == conds]
    
    
    mids_less_zero <- sum(spontaneous_histo$mids < 0)
    # if it is significant determine onset 
    if(signif_table$signif[signif_table$condition == conds]){
      
      # CPM changepoints
      change_points <- processStream(x=full_histograms[[conds]], cpmType = "Mann-Whitney",
                                     ARL0 = 5000, startup = mids_less_zero)
      if(length(change_points$changePoints)!= 0){
        chnpoints <- change_points$changePoints[change_points$changePoints > mids_less_zero]
        if(length(chnpoints)!=0){
          times <- spontaneous_histo$mids[chnpoints]
          cell_changepoints <- rbind(cell_changepoints, data.frame(cell=rep(neuron, length(times)), 
                                                                   condition=rep(conds, length(times)), 
                                                                   changepoints=times))
          
      }
    }
  }
  
  # does the cell show the usual golgi depression
  spont_analysis$depression_to_stim[spont_analysis$cell == neuron] <- any(unlist(response_analysis[response_analysis$cell==neuron,
                                                                                                   c("forelimb", "hindlimb", "vib")]))
  
  plot_isi_df <- filter(cell_isis, cell==neuron)
  #kruskal wallis rank sum on cell isi for this cell, save in text file
  kruskal_test <- kruskal.test(isis~condition, data = plot_isi_df)
  sink(paste0(neuron, " kruskal wallis test of isis.txt"))
  print(kruskal_test)
  sink()
  
  mcp5_response <- response_analysis$MCP1[response_analysis$cell==neuron]
  mcp1_response <- response_analysis$MCP5[response_analysis$cell==neuron]
  # any response to MCP stims
  if((!is.na(mcp5_response) & mcp5_response) |(!is.na(mcp1_response) & mcp1_response)){
    kruskal_isi$response_to_mcp[kruskal_isi$cell == neuron] <- T
    cell_isis$signif[cell_isis$cell == neuron] <- T
    plot_isi_df$signif <- T
  } else {
    kruskal_isi$response_to_mcp[kruskal_isi$cell == neuron] <- F
    cell_isis$signif[cell_isis$cell == neuron] <- F
    plot_isi_df$signif <- F
  }
  # store p value
  kruskal_isi$p_value[kruskal_isi$cell == neuron] <- kruskal_test$p.value
  
  # plots of raw isi
  # find limits of plot
  lower_whisker <- numeric(0)
  upper_whisker <- numeric(0)
  for(conditions in unique(plot_isi_df$condition)){
    df <- filter(plot_isi_df, condition==conditions)
    box_stats <- boxplot.stats(df$isis)
    lower_whisker <- c(lower_whisker, box_stats$stats[1])
    upper_whisker <- c(upper_whisker, box_stats$stats[5])
  }
  ylim1 <- c(min(lower_whisker), max(upper_whisker)) 
  
  # rotate labels, change y values to ms
  kruskal_plot <- ggplot(plot_isi_df, aes(x=condition, y=isis, fill=signif)) + geom_boxplot(outlier.colour = NA) +
    scale_y_continuous(limits=ylim1*1.05, labels = function(x)x*1000, name = "interspike interval (ms)") +
    ggtitle(ifelse(kruskal_test$p.value < 0.05, "p < 0.05", "p > 0.05")) + guides(fill=guide_legend(title="Response to MCP Stimulation?")) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 12))
  
  # save 
  ggsave(kruskal_plot, filename = paste0(neuron, " plots of isis for spontaneous and after MCP.pdf"))
  
  # if there is signif 1MCP, what are the time of driven spikes (if it is)
  if((!is.na(mcp1_response) & mcp1_response)){
    spike_table <- read.table(grep("1MCP", files, value = T), header = T, sep = "\t",
                              colClasses = c("numeric", "numeric", "character"), fill = T)
    spike_times <- as.numeric(unlist(strsplit(spike_table[, 3], split = ",")))
    
    # spike times after stim
    spike_times <- spike_times[spike_times > 0 & spike_times < 0.005]
    spike_plot <- ggplot() + geom_histogram(aes(x=spike_times), binwidth = 0.0001) + 
      theme_classic() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(limits=c(0, 0.005),labels=function(x)x*1000,
                                                                               name = "latency of spikes in response\nto MCP stimulus (ms)",
                                                                               expand = c(0,0))
    ggsave(spike_plot, filename = paste0(neuron, " spike times in response to 1MCP.pdf"))

    # rbind spike times
    driven_spikes <- rbind(driven_spikes, data.frame(cell=rep(neuron, length(spike_times)), times=spike_times))
  }
  
  # if response to 5MCP, work out probabilty of spike to first, second, third...
  if((!is.na(mcp5_response) & mcp5_response)){
    spike_table <- read.table(grep("5MCP", files, value = T), header = T, sep = "\t",
                              colClasses = c("numeric", "numeric", "character"), fill = T)
    n_sweeps <- length(spike_table[, 3])
    spike_times <- as.numeric(unlist(strsplit(spike_table[, 3], split = ",")))
    
    spikes_to_stim <- c(length(spike_times[spike_times > 0 & spike_times < 0.0033])/n_sweeps,
                        length(spike_times[spike_times > 0.0033 & spike_times < 0.0066])/n_sweeps,
                        length(spike_times[spike_times > 0.0066 & spike_times < 0.0099])/n_sweeps,
                        length(spike_times[spike_times > 0.0099 & spike_times < 0.0132])/n_sweeps,
                        length(spike_times[spike_times > 0.0132 & spike_times < 0.0165])/n_sweeps)
    
    driven_probability <- rbind(driven_probability, data.frame(cell=rep(neuron, 5), stim_number = seq(1,5,1),
                                                               probability=spikes_to_stim))
    
  }
  setwd("..")
}

