# Analysis of pons cells
# 1> spontaneous analysis: Hz, ISIs, range ISI
# 2> limb stimulation
# 3> antidromic inhibition

# libraries needed
library(ggplot2)
library(dplyr)
library(reader)
library(cpm)

# functions
source("functions/multivariate_gaussian.R")
source("functions/histogram_func.R")

# standard error function 
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# parameters for histograms
pre <- -0.2
post <- 0.4
binwidth <- 0.01

# choose the correct directory
setwd(choose.dir())

# # list all folders
folders <- list.dirs()
folders <- gsub("./", x = folders, replacement = "", fixed = T)
folders <- folders[folders!="."]
cell_num <- length(folders)

# spontaneous analysis
spont_analysis <- data.frame(cell = folders, Hz=rep(NA, cell_num), mean_ISI = rep(NA, cell_num),
                             med_ISI = rep(NA, cell_num), min_ISI = rep(NA, cell_num), max_ISI = rep(NA, cell_num),
                             cv_ISI = rep(NA, cell_num))

# limb stimulation results
response_analysis <- data.frame(cell = folders, contra_forelimb = rep(NA, cell_num), ipsi_forelimb = rep(NA, cell_num),
                           contra_hindlimb = rep(NA, cell_num), ipsi_hindlimb=rep(NA, cell_num))

# onset data frame
onset_df <- data.frame(cell= folders, contra_forelimb = rep(NA, cell_num), ipsi_forelimb = rep(NA, cell_num),
                       contra_hindlimb = rep(NA, cell_num), ipsi_hindlimb=rep(NA, cell_num))

# offset data frame
offset_df <- data.frame(cell = folders , contra_forelimb = rep(NA, cell_num), ipsi_forelimb = rep(NA, cell_num),
                        contra_hindlimb = rep(NA, cell_num), ipsi_hindlimb=rep(NA, cell_num))

# antidromic probability to stim
antidromic_df  <- data.frame(cell = rep(folders, each =3), stim_number=rep(1:3, cell_num),
                             probability = rep(NA, cell_num*3))

# cell isis
cell_isis <- as.data.frame(matrix(nrow = 0, ncol=2))
colnames(cell_isis) <- c("cell", "isi")

# timing to antidromic
antidromic_timing <- as.data.frame(matrix(nrow=0, ncol=3))
colnames(antidromic_timing) <- c("cell", "stim_number", "latency")

# cell changepoints
cell_changepoints <- as.data.frame(matrix(nrow=0, ncol=3))
colnames(cell_changepoints) <- c("cell", "conditions", "changepoints")

# colours for stimulation
colour_key <- list('contra forelimb'="#ff0000", 'contra hindlimb'="#1d0eba", 'ipsi forelimb'="#8a10de", 'ipsi hindlimb'="#f4a300")

# correction table, some antidromic have delays and interpulse is varied
correction_table <- data.frame(cell=c("cell 13", "cell 18", "cell 19",
                                      "cell 20", "cell 21", "cell 22",
                                      "cell 29", "cell 30", "cell 31",
                                      "cell 35", "cell 38", "cell 45",
                                      "cell 46", "cell 49", "cell 51",
                                      "cell 53", "cell 56"),
                               delay = c(0.02, 0, 0, 0.002, 0.0003,
                                         0, 0.003, 0.003, 0.0014,
                                         0.0012, 0.0022, 0.0022, 0,
                                         0, 0, 0, 0.003),
                               interpulse = c(rep(0.0043, 8), 0.0038,
                                              rep(0.0045, 8)))

cell_sweeps <- data.frame(cell = c("cell 31", "cell 32", "cell 33",
                                   "cell 34", "cell 35", "cell 36",
                                   "cell 37", "cell 38", "cell 40",
                                   "cell 41", "cell 42", "cell 43",
                                   "cell 44", "cell 45", "cell 46",
                                   "cell 47", "cell 48", "cell 49",
                                   "cell 50", "cell 51", "cell 52",
                                   "cell 53", "cell 54", "cell 55",
                                   "cell 56"), 
                          n_sweep = c(140, 160, 570, 320,
                                      500, 340, 540, 520,
                                      440, 360, 530, 520,
                                      730, 520, 590, 500,
                                      520, 520, 540, 520,
                                      650, 560, 680, 590,
                                      560)) 

# iterate through folders
for(neuron in folders){
  # move in to cell directory
  setwd(neuron)
  files <- list.files()
  cat("working on", which(folders==neuron), "out of", length(folders), "\n")
  spontaneous_file <- paste0(neuron, ".txt")
  
  # every cell has spontaneous, analyse this first
  total_file_time <- as.numeric(readLines(grep(spontaneous_file, files, value = T))[5])
  spont_times <- as.numeric(readLines(grep(spontaneous_file, files, value = T))[6:file.nrow(grep(spontaneous_file, files, value = T))])
  
  # stats for spont; Hz
  spont_analysis$Hz[spont_analysis$cell == neuron] <- length(spont_times)/total_file_time
  
  # stats for spont; ISIs, CV log ISI
  shifted_spont <- c(spont_times[2:length(spont_times)], NA)
  ISIs <- na.omit(shifted_spont - spont_times)
  
  # ISI stats, mean, median, min, max, cv
  spont_analysis$mean_ISI[spont_analysis$cell == neuron] <- mean(ISIs)
  spont_analysis$med_ISI[spont_analysis$cell == neuron] <- median(ISIs)
  spont_analysis$min_ISI[spont_analysis$cell == neuron] <- min(ISIs)
  spont_analysis$max_ISI[spont_analysis$cell == neuron] <- max(ISIs)
  spont_analysis$cv_ISI[spont_analysis$cell == neuron] <- sd(ISIs)/mean(ISIs)
  
  # plot of ISIs
  spontISI <- ggplot() + geom_histogram(aes(x=ISIs), breaks =seq(0, max(ISIs), 0.01)) + 
    theme_classic() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(labels=function(x)x*1000,
                                                                             name = "interspike interval (ms)",
                                                                             expand = c(0,0))
  # create a directory for the ISI plot
  dir.create("spontaneous isi plot")
  setwd("spontaneous isi plot")
  ggsave(spontISI, filename = "spontaneous ISI long.pdf")
  shortISI <- spontISI + coord_cartesian(xlim=c(0, 1))
  ggsave(shortISI, filename = "spontaneous ISI short.pdf")
  setwd("..")
  
  # bind the ISIs to the cell isi data frame
  cell_isis <- rbind(cell_isis, data.frame(cell = rep(neuron, length(ISIs)), 
                                           isi = ISIs))
  
  # see if limb stim files exist
  stim_conditions <- c("contra forelimb", "contra hindlimb", "ipsi forelimb", "ipsi hindlimb")
  greps <- logical()
  for(stims in stim_conditions){
    greps <- c(greps, any(grepl(stims, x = files)))
  }
  
  # if limb stimulation exist
  if(any(greps)){
    # as true create stimulus histogram
    spontaneous_histo <- spontaneousHist(spont_times = spont_times, pre = pre, post = post, iteration = 500,
                                         sweeps=cell_sweeps$n_sweep[cell_sweeps$cell == neuron])
    spontaneous_histo$full <- na.omit(spontaneous_histo$full)
    spontaneous_histo$post <- na.omit(spontaneous_histo$post)
    smooth_spont <- as.data.frame(matrix(nrow=length(spontaneous_histo$mids)*nrow(spontaneous_histo$full), ncol=6))
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
    
    # conditions that have files
    true_stims <- which(greps==T)
    conditions <- stim_conditions[true_stims]
    for(cond in conditions){
      
      # read in the spike times
      spike_table <- read.table(grep(cond, files, value = T), header = T, sep = "\t",
                                colClasses = c("numeric", "numeric", "character"), fill = T)
      spike_times <- as.numeric(unlist(strsplit(spike_table[, 3], split = ",")))
      
      # process the spike times
      stim_hist <- stimulusHistogram(spike_times = spike_times, pre = pre, post = post)
      
      # place it in the stimulus list
      stimulus_histograms[[cond]] <- stim_hist$post
      full_histograms[[cond]] <- stim_hist$full
      
      value <- smooth.spline(x=spontaneous_histo$mids, y=stim_hist$full, spar=0.3)$y
      value2 <- smooth.spline(x=spontaneous_histo$mids, y=stim_hist$full, spar=0.5)$y
      
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
      
      # stimulation plot
      stim_plot <- ggplot(smooth_stim, aes(x=time, y=value2, group=group, colour=condition)) +
        geom_line(aes(alpha=alpha), size=1) + scale_alpha_identity() + theme_classic() +
        scale_x_continuous(name="time relative to stimulus (ms)", limits=c(pre, post), labels=function(x)x*1000,
                           expand = c(0,0)) + 
        scale_y_continuous(name="% activity compared to baseline", limits=c(0,y_max), expand=c(0,0)) + 
        scale_color_manual(values = c(colour_key[[cond]], "#000000"), guide=guide_legend(title = ""), labels=c(cond, "spontaneous"))
      
      dir.create("stimulation plots")
      setwd("stimulation plots")
      ggsave(stim_plot, filename = paste0(cond, " stimulation plot.pdf"))
      setwd("..")
    }
    
    # compare the x1..xn bins to spontaneous x1 ... xn bins
    params <- normalParams(spontaneous_histo$post)
    signif_table <- gaussianConfInt(x = stimulus_histograms, mu = params$mu,
                                    Sigma = params$sigma, alpha = 0.0001)
    
    # iterate through the significance table
    for(conds in signif_table$condition){
      con <- paste(strsplit(conds, " ")[[1]], collapse = "_")
      
      # place the value in the response table
      response_analysis[response_analysis$cell==neuron, con] <- signif_table$signif[signif_table$condition == conds]
      
      # mids less than zero, number needed as the startup for the changepoint analysis
      mids_less_zero <- sum(spontaneous_histo$mids < 0)
      
      if(signif_table$signif[signif_table$condition == conds]){
        
        # CPM changepoints
        change_points <- processStream(x=full_histograms[[conds]], cpmType = "Mann-Whitney",
                                       ARL0 = 1000, startup = mids_less_zero)
        if(length(change_points$changePoints)!= 0){
          chnpoints <- change_points$changePoints[change_points$changePoints > mids_less_zero]
          if(length(chnpoints)!=0){
            times <- spontaneous_histo$mids[chnpoints]
            cell_changepoints <- rbind(cell_changepoints, data.frame(cell=rep(neuron, length(times)), 
                                                                     condition=rep(conds, length(times)), 
                                                                     changepoints=times))
            onset_df[onset_df$cell == neuron, con] <- times[1]
            
          }
        }
      }
    }
  }
  
  # MCP analysis 
  if(any(grepl("MCP", x = files))){
    # read in the spike times
    spike_table <- read.table(grep("MCP", files, value = T), header = T, sep = "\t",
                              colClasses = c("numeric", "numeric", "character"), fill = T)
    spike_times <- as.numeric(unlist(strsplit(spike_table[, 3], split = ",")))
    n_sweeps <- length(spike_table[, 3])
    
    # correct the delay and define the interpulse interval
    delay <- correction_table$delay[correction_table$cell == neuron]
    interpulse <- correction_table$interpulse[correction_table$cell == neuron]
    
    # correct the spike times
    spike_times <- spike_times - delay
    
    # split into first, second and third
    first_stim <- spike_times[spike_times > 0 & spike_times <= interpulse]
    second_stim <- spike_times[spike_times > interpulse & spike_times <= interpulse*2]
    third_stim <- spike_times[spike_times > interpulse*2 & spike_times <= interpulse*3]
    
    # latency from pulses
    first_latency <- first_stim
    second_latency <- second_stim - interpulse
    third_latency <- third_stim - interpulse*2
    
    # probability to pulses
    first_probability <- length(first_stim)/n_sweeps
    second_probability <- length(second_stim)/n_sweeps
    third_probability <- length(third_stim)/n_sweeps
    
    # add values to antidromic dataframes
    antidromic_df$probability[antidromic_df$cell == neuron & 
                                    antidromic_df$stim_number==1] <- first_probability
    antidromic_df$probability[antidromic_df$cell == neuron & 
                                antidromic_df$stim_number==2] <- second_probability
    antidromic_df$probability[antidromic_df$cell == neuron & 
                                antidromic_df$stim_number==3] <- third_probability
    
    # and the antidromic latency
    antidromic_timing <- rbind(antidromic_timing, data.frame(cell=rep(neuron, length(first_latency) +
                                                                        length(second_latency) + length(third_latency)),
                                                             stim_number=c(rep(1, length(first_latency)),
                                                                           rep(2, length(second_latency)),
                                                                           rep(3, length(third_latency))),
                                                             latency = c(first_latency, second_latency, third_latency)))
  } 
  # back to main directory
  setwd("..")
}
write.csv(spont_analysis, "spontaneous analysis.csv", row.names = F)
write.csv(response_analysis, "response analysis.csv", row.names = F)
write.csv(onset_df, "onset dataframe.csv", row.names = F)
write.csv(offset_df, "offset dataframe.csv", row.names = F)
write.csv(antidromic_df, "antidromic probability dataframe.csv", row.names = F)
write.csv(antidromic_timing, "antidromic timing.csv", row.names = F)
write.csv(cell_isis,"cell isi.csv", row.names = F)
write.csv(cell_changepoints, "changepoints.csv", row.names = F)


# spontaneous analysis
Hz_hist <- ggplot(spont_analysis, aes(x=Hz)) + geom_histogram(breaks=seq(0, max(spont_analysis$Hz)+1, 1)) +
  theme_classic() + scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0))
ggsave(Hz_hist, filename = "spontanoeus Hz histogram.pdf")


# comparision of average ISI; means and medians
med_mean <- data.frame(average=c(rep("Mean ISI",nrow(spont_analysis)), rep("Median ISI", nrow(spont_analysis))),
                       values = c(spont_analysis$mean_ISI, spont_analysis$med_ISI))

md_mn_hist <- ggplot(med_mean, aes(x=values)) + geom_histogram(breaks=seq(0, 4.5, 0.1)) + facet_wrap(~average, ncol=1) +
  xlab("Interspike Interval (s)")
ggsave(md_mn_hist, filename = "comparison of median and mean ISI.pdf")

# min isi
min_isi_plot <- ggplot(spont_analysis, aes(x=min_ISI)) + geom_histogram(breaks=seq(0, max(spont_analysis$min_ISI)+0.01, 0.001)) +
  theme_classic() + scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
  xlab("Minimum Interspike Interval (s)")
ggsave(min_isi_plot, filename = "minimum isi histogram.pdf")

max_isi_plot <- ggplot(spont_analysis, aes(x=max_ISI)) + geom_histogram(breaks=seq(0, max(spont_analysis$max_ISI)+1, .5)) +
  theme_classic() + scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
  xlab("Maximum Interspike Interval (s)")
ggsave(max_isi_plot, filename = "maximum isis histogram.pdf")


# response analysis, remove complete absences
row_remove <- numeric()
for(i in 1:nrow(response_analysis)){
  if(sum(is.na(response_analysis[i, ])) == 4){
    row_remove <- c(row_remove, i)
  }
}

# response table with only the tested cells
tested_response <- response_analysis[-(row_remove), ]

# deal with the NAs by making them false
t_response <- tested_response
t_response[is.na(tested_response)] <- F

conditions <- character()
for(i in 1:nrow(t_response)){
  conditions <- c(conditions, paste(t_response[i, 2:5], collapse = " "))
}

# used table conditions to count the unique incidents 
response_count <- data.frame(Responses=c("Contralateral Forelimb", "Contra & Ipsilateral Forelimb", "Contra & Ipsilateral Hindlimb",
                                         "Contra & Ipsilateral Forelimb + Ipsilateral Hindlimb", "All Limbs", "No Response", "TOTAL"),
                             Counts=c(3, 9, 1, 3, 2, 7, 25))

library(grid)
library(gridExtra)
ggsave(tableGrob(response_count, rows = NULL), filename="count table of limb responses.pdf")

# where a limb has been tested what proportion have a response
co_fl <- signif(sum(na.omit(tested_response$contra_forelimb))/length(na.omit(tested_response$contra_forelimb)) * 100, 4)
ip_fl <- signif(sum(na.omit(tested_response$ipsi_forelimb))/length(na.omit(tested_response$ipsi_forelimb)) * 100, 4)
co_hl <- signif(sum(na.omit(tested_response$contra_hindlimb))/length(na.omit(tested_response$contra_hindlimb)) * 100, 4)
ip_hl <- signif(sum(na.omit(tested_response$ipsi_hindlimb))/length(na.omit(tested_response$ipsi_hindlimb)) * 100, 4)

response_percentages <- data.frame('Stimulus condition'=c("Contralateral Forelimb", "Contralateral Hindlimb",
                                                          "Ipsilateral Forelimb", "Contralateral Hindlimb"),
                                   'Percent of Responses'=c(co_fl, co_hl, ip_fl, ip_hl))
ggsave(tableGrob(response_percentages, rows = NULL), filename = "percentage of responses to conditions.pdf")


# antidromic probabilities and timing
anti_prob <- na.omit(antidromic_df)
anti_prob_plot <- ggplot(anti_prob, aes(x=stim_number, y=probability, colour=cell)) + geom_point() + geom_line() + 
  scale_x_continuous(breaks=c(1, 2, 3), name = "Antidromic MCP Pulse Number") + ylab("Probability of an Antidromic Spike")
ggsave(anti_prob_plot, filename = "probability of antidromic spike.pdf")


# anti spikes
anti_spikes <- read.csv(file.choose())
antidromic_spike <- ggplot(anti_spikes, aes(x=delay, y=spike_number, group=cell)) + geom_line() + geom_point() + 
  scale_x_discrete(limits=c("No Delay", "Delay"), name="Delay from spontaneous spike\nto antidromic stimulus onset") +
  ylab("Number of antidromic spikes")
ggsave(antidromic_spike, filename = "antidromic spike counts.pdf")

# timing of antidromic spikes
anti_timing <- ggplot(antidromic_timing, aes(x=latency)) + geom_histogram(breaks=seq(0, 0.004, 0.0001)) + facet_wrap(~cell, ncol = 3, scales = "free_y") + 
  scale_x_continuous(labels=function(x)x*1000, name="Antidromic Spike Latency (ms)") + ylab("Count")
ggsave(anti_timing, filename = "antidromic spike latencies.pdf")

avg_timings <- antidromic_timing %>%
  group_by(cell) %>%
  summarise(avg=mean(latency))
avg_boxplot <- ggplot(avg_timings, aes(x="Average\nLatencies" ,y=avg)) + geom_boxplot() + coord_flip() +
  xlab("") + scale_y_continuous(labels=function(x)x*1000, name="Average Antidromic Spike Latencies (ms)", limits = c(0, NA))
ggsave(avg_boxplot, filename = "average antidromic latencies boxplot.pdf")
