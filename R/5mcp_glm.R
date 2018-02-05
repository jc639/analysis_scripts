# 5mcp spike probability glm
# driven spike = [0 | 1]
# predictors = previous spike driven, subsequent spike driven, stim number
library(ggplot2)
library(dplyr)
response_analysis <- read.csv(file.choose())

# which ones are significant
signif_mcp5 <- response_analysis$cell[response_analysis$MCP5==T & !is.na(response_analysis$MCP5)]

# set the wd to the text file folder
setwd(choose.dir())

# create empty dataframe for glm
driven_prob <- as.data.frame(matrix(nrow = 0, ncol = 5)) 
colnames(driven_prob) <- c("cell", "stim_number", "previous_spike", "subsequent_spike", "spike")

# iterate through cells
for(cell in signif_mcp5){
  setwd(as.character(cell))
  
  # list files
  files <- list.files()
  
  # read in the spike table
  spike_table <- read.table(grep("5MCP", files, value = T), header = T, sep = "\t",
                            colClasses = c("numeric", "numeric", "character"), fill = T)
  
  # split out into sweeps
  spike_sweeps <- lapply(strsplit(spike_table[, 3], split=","), as.numeric)
  
  # iterate through sweeps
  for(i in 1:length(spike_sweeps)){
    # spike times from sweep
    spike_times <- spike_sweeps[[i]]
    
    # length of any spike times between stimulus pulses
    first_spike <- length(spike_times[spike_times > 0 & spike_times < 0.0033])
    second_spike <- length(spike_times[spike_times > 0.0033 & spike_times < 0.0066])
    third_spike <- length(spike_times[spike_times > 0.0066 & spike_times < 0.0099])
    fourth_spike <- length(spike_times[spike_times > 0.0099 & spike_times < 0.0132])
    fifth_spike <- length(spike_times[spike_times > 0.0132 & spike_times < 0.0165])
    
    # if first spikes exist, add a row to driven_prob with spike indicated as 1, else not
    # indicate a spike on previous or subsequent with a y/n
    if(first_spike > 0){
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 1, previous_spike = NA, subsequent_spike = ifelse(second_spike > 0, "yes", "no"),
                                    spike = 1))
    } else {
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 1, previous_spike = NA, subsequent_spike = ifelse(second_spike > 0, "yes", "no"),
                                                   spike = 0))
    }
    
    if(second_spike > 0){
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 2, previous_spike = ifelse(first_spike > 0, "yes", "no"), subsequent_spike = ifelse(third_spike > 0, "yes", "no"),
                                                   spike = 1))
    } else {
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 2, previous_spike = ifelse(first_spike > 0, "yes", "no"), subsequent_spike = ifelse(third_spike > 0, "yes", "no"),
                                                   spike = 0))
    }
    
    if(third_spike > 0){
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 3, previous_spike = ifelse(second_spike > 0, "yes", "no"), subsequent_spike = ifelse(fourth_spike > 0, "yes", "no"),
                                                   spike = 1))
    } else {
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 3, previous_spike = ifelse(second_spike > 0, "yes", "no"), subsequent_spike = ifelse(fourth_spike > 0, "yes", "no"),
                                                   spike = 0))
    }
    
    if(fourth_spike > 0){
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 4, previous_spike = ifelse(third_spike > 0, "yes", "no"), subsequent_spike = ifelse(fifth_spike > 0, "yes", "no"),
                                                   spike = 1))
    } else {
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 4, previous_spike = ifelse(third_spike > 0, "yes", "no"), subsequent_spike = ifelse(fifth_spike > 0, "yes", "no"),
                                                   spike = 0))
    }
    
    if(fifth_spike){
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 5, previous_spike = ifelse(fourth_spike > 0, "yes", "no"), subsequent_spike = NA,
                                                   spike = 1))
    } else {
      driven_prob <- rbind(driven_prob, data.frame(cell = cell, stim_number = 5, previous_spike = ifelse(fourth_spike > 0, "yes", "no"), subsequent_spike = NA,
                                                   spike = 0))
    }
  }
  setwd("..")
}
# glms of probabilites, two models one for stims 1-5; one for 2-4, does previous or subsequent spike effect probability
stim_num_glm <- glm(spike~as.factor(stim_number), data=driven_prob, family = binomial(link = "logit"))
sink("glm of 5mcp stim number probability.txt")
summary(stim_num_glm)
sink()

# create dataframe of probability
stim_num_df <- data.frame(stim_num = 1:5, probability = predict(stim_num_glm, data.frame(stim_number=as.factor(1:5)), 
                                                                type = "response"))

# individual cell probabilities
cell_stim_prob <- driven_prob %>%
  group_by(cell, stim_number) %>%
  summarise(probability = sum(spike)/length(spike))

# plot of probability by spike number
stim_number_plot <- ggplot() + geom_bar(data = stim_num_df, aes(x=stim_num, y=probability), stat = "identity") +
  geom_point(data=cell_stim_prob, aes(x=stim_number, y=probability), colour="blue", alpha=0.3) + 
  ylab("Probability of Spike\nFollowing MCP Stimulus Pulse") + xlab("MCP Stimulus Pulse Number")
ggsave(stim_number_plot, filename = "probability of spike following stimulus pulse.pdf")


##### same again but without cell 22 - outlier ####################
# glms of probabilites, two models one for stims 1-5; one for 2-4, does previous or subsequent spike effect probability
stim_num_glm <- glm(spike~as.factor(stim_number), data=driven_prob, subset=cell!="cell 22", family = binomial(link = "logit"))
sink("glm of 5mcp stim number probability minus cell 22.txt")
summary(stim_num_glm)
sink()

# create dataframe of probability
stim_num_df <- data.frame(stim_num = 1:5, probability = predict(stim_num_glm, data.frame(stim_number=as.factor(1:5)), 
                                                                type = "response"))

# individual cell probabilities
cell_stim_prob <- driven_prob %>%
  filter(cell!="cell 22") %>%
  group_by(cell, stim_number) %>%
  summarise(probability = sum(spike)/length(spike))

# plot of probability by spike number
stim_number_plot <- ggplot() + geom_bar(data = stim_num_df, aes(x=stim_num, y=probability), stat = "identity") +
  geom_point(data=cell_stim_prob, aes(x=stim_number, y=probability), colour="blue", alpha=0.3) + 
  ylab("Probability of Spike\nFollowing MCP Stimulus Pulse") + xlab("MCP Stimulus Pulse Number")
ggsave(stim_number_plot, filename = "probability of spike following stimulus pulse minus cell 22.pdf")

##########################################################################################


# probability if previous pulse has generated a spike
filt_driven_prob <- filter(driven_prob, stim_number%in%2:4)

# whether neighbouring stim has generated a spike causing a spike
neigh_glm <- glm(spike~as.factor(stim_number)*previous_spike*subsequent_spike, data=filt_driven_prob,
                 family = binomial(link = "logit"))


neigh_prob <- data.frame(stim_number=as.factor(c(rep(2,4), rep(3, 4), rep(4, 4))),
                         previous_spike = rep(c("no", "no", "yes", "yes"), 3),
                         subsequent_spike = rep(c("no", "yes", "no", "yes"), 3),
                         probability=rep(NA, 12))

  
neigh_prob$probability <- predict(neigh_glm, neigh_prob[, -4], type = "response")
neigh_prob$condition <- apply(neigh_prob[,c(2,3)], MARGIN = 1, paste, collapse=" - ")

filt_driven_prob$condition <- apply(filt_driven_prob[, c(3,4)], 1, paste, collapse=" - ")

cond_prob <- filt_driven_prob %>%
  group_by(cell, stim_number, condition) %>%
  summarise(probability=sum(spike)/length(spike))

stim_number_plot <- ggplot() + geom_bar(data=neigh_prob, aes(x=condition, y=probability), stat = "identity") +
  geom_point(data=cond_prob, aes(x=condition, y=probability), colour="blue", alpha = 0.3) + facet_wrap(~stim_number) +
  xlab("Spike to Previous Stimulus Pulse? -- Spike to Subsequent Stimulus Pulse?") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Probability of Response to Stimulus") + ggtitle("Stimulus Pulse Number")

sink("glm of spike probability to stim number.txt")
summary(neigh_glm)
sink()

ggsave(stim_number_plot, filename = "spike probability.pdf")


#### same but without cell 22 - outlier #############


# probability if previous pulse has generated a spike
filt_driven_prob <- filter(driven_prob, stim_number%in%2:4, cell!="cell 22")

# whether neighbouring stim has generated a spike causing a spike
neigh_glm <- glm(spike~as.factor(stim_number)*previous_spike*subsequent_spike, data=filt_driven_prob,
                 family = binomial(link = "logit"))


neigh_prob <- data.frame(stim_number=as.factor(c(rep(2,4), rep(3, 4), rep(4, 4))),
                         previous_spike = rep(c("no", "no", "yes", "yes"), 3),
                         subsequent_spike = rep(c("no", "yes", "no", "yes"), 3),
                         probability=rep(NA, 12))


neigh_prob$probability <- predict(neigh_glm, neigh_prob[, -4], type = "response")
neigh_prob$condition <- apply(neigh_prob[,c(2,3)], MARGIN = 1, paste, collapse=" - ")

filt_driven_prob$condition <- apply(filt_driven_prob[, c(3,4)], 1, paste, collapse=" - ")

cond_prob <- filt_driven_prob %>%
  group_by(cell, stim_number, condition) %>%
  summarise(probability=sum(spike)/length(spike))

stim_number_plot <- ggplot() + geom_bar(data=neigh_prob, aes(x=condition, y=probability), stat = "identity") +
  geom_point(data=cond_prob, aes(x=condition, y=probability), colour="blue", alpha = 0.3) + facet_wrap(~stim_number) +
  xlab("Spike to Previous Stimulus Pulse? -- Spike to Subsequent Stimulus Pulse?") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Probability of Response to Stimulus") + ggtitle("Stimulus Pulse Number")

sink("glm of spike probability to stim number without cell 22.txt")
summary(neigh_glm)
sink()

ggsave(stim_number_plot, filename = "spike probability without cell 22.pdf")
