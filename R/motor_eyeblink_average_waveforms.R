# motor eyeblink waveform average:
#   1 - average of every 20 trials for individuals from 1 - 120
#   2 - average of the different button types for each id over 120 trials
#   3 - average of before switching hands and after
#   4 - average of before experimenter presses and after

# libraries
library(dplyr)
library(ggplot2)

# se function
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# read in trial waveform csv
setwd(choose.dir())

# trial waveform
trial_wv <- read.csv("trial waveform.csv")

# create folder for plots
dir.create("average plots")
setwd("average plots")

# change block names
trial_wv$blocks <- gsub("(", x=trial_wv$blocks, replacement = "", fixed=T)
trial_wv$blocks <- gsub("[", x=trial_wv$blocks, replacement = "", fixed=T)
trial_wv$blocks <- gsub("]", x=trial_wv$blocks, replacement = "", fixed=T)
trial_wv$blocks <- gsub(",", x=trial_wv$blocks, replacement = "-", fixed=T)

# correct blocks
training_blocks <- unique(trial_wv$blocks)[1:6]
switch_blocks <- unique(trial_wv$blocks)[6:7]
outside_blocks <- unique(trial_wv$blocks)[10:11]

# filter to just one group, create directory
for(grp in c("control", "experimental")){
  group_df <- trial_wv %>%
    filter(group==grp, response%in%c("CR", "No_CR"))
  
  # create group directory and move in to it
  dir.create(grp)
  setwd(grp)
  
  # unique ids of subjects
  uni_id <- unique(group_df$id)
  
  # create folder
  for(ident in uni_id){
    dir.create(as.character(ident))
    setwd(as.character(ident))
    
    id_df <- group_df %>%
      filter(id==ident)
    
    # data frame for the first training period
    blocks_df <- id_df %>%
      filter(blocks%in%training_blocks) %>%
      group_by(blocks, time) %>%
      summarise(avg = mean(value)*-1,
                se = se(value),
                alp = ifelse(unique(blocks)%in%c("1-20","100-120"), 1, 0.7))
    
    # relevel
    blocks_df$blocks <- factor(as.character(blocks_df$blocks), levels = training_blocks)
    
    # grouped plot all together
    grouped_plot <- ggplot(blocks_df, aes(x=time, y=avg, group=blocks)) + geom_line(aes(colour=blocks), size=1) +
      theme_classic() + geom_vline(xintercept = 0, linetype=2, alpha=0.4) + scale_alpha_identity() +
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks = c(-0.5,0,0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.05) + ylab("blink response (arbitary V)") +
      scale_color_brewer(palette = "RdYlBu")
    
    # facet plot of the 20 trial blocks
    facet_plot <- ggplot(blocks_df, aes(x=time, y=avg, group=blocks)) + geom_line(aes(colour=blocks), size=1) +
      theme_classic() + geom_vline(xintercept = 0, linetype=2, alpha=0.4) +
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks = c(-0.5,0,0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.3) + ylab("blink response (arbitary V)") +
      facet_wrap(~blocks, strip.position = "left")
    
    # button data frame for the training period
    button_df <- id_df %>%
      filter(blocks%in%training_blocks[2:6]) %>%
      group_by(button, time) %>%
      summarise(avg = mean(value)*-1,
                se = se(value))
    
    # renane the test periods
    button_df$button <- gsub(pattern = "BTest", x = button_df$button, replacement = "B test")
    button_df$button <- gsub(pattern = "CTest", x = button_df$button, replacement = "C test")
    
    # button plot, faceted
    button_plot <- ggplot(button_df, aes(x=time, y=avg, group=button)) + geom_line(aes(colour=button), size=1) +
      theme_classic() + geom_vline(xintercept=0, linetype=2, alpha=0.4) + 
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks=c(-0.5, 0, 0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.3) + ylab("blink response (arbitary V)") +
      facet_wrap(~button)
    
    # switch hands dataframe
    switch_df <- id_df %>%
      filter(blocks%in%switch_blocks) %>%
      group_by(blocks, time) %>%
      summarise(avg = mean(value)*-1,
                se = se(value),
                alp=ifelse(unique(blocks)=="120-140", 1, 0.8))
    
    # switch plot grouped together
    switch_plot <- ggplot(switch_df, aes(x=time, y=avg, group=blocks)) + geom_line(aes(colour=blocks, alpha=alp), size=1) +
      theme_classic() + geom_vline(xintercept=0, linetype=2, alpha=0.4) + scale_alpha_identity() +
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks=c(-0.5, 0, 0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.1) + ylab("blink response (arbitary V)") +
      scale_colour_discrete(labels=c("right hand", "left hand"), name="switching hands")
    
    # switch plot faceted
    switch_facet <- ggplot(switch_df, aes(x=time, y=avg, group=blocks)) + geom_line(aes(colour=blocks), size=1) +
      theme_classic() + geom_vline(xintercept=0, linetype=2, alpha=0.4) + 
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks=c(-0.5, 0, 0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.3) + ylab("blink response (arbitary V)") +
      scale_colour_discrete(labels=c("right hand", "left hand"), name="switching hands") + facet_wrap(~blocks, ncol=2, 
                                                                                                      strip.position = "left")
    # outsider pressing dataframe
    outside_df <- id_df %>%
      filter(blocks%in%outside_blocks) %>%
      group_by(blocks, time) %>%
      summarise(avg = mean(value)*-1,
                se = se(value),
                alp=ifelse(unique(blocks)=="200-220", 1, 0.8))
    
    # outside plotm together
    outside_plot <- ggplot(outside_df, aes(x=time, y=avg, group=blocks)) + geom_line(aes(colour=blocks, alpha=alp), size=1) +
      theme_classic() + geom_vline(xintercept=0, linetype=2, alpha=0.4) + scale_alpha_identity() +
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks=c(-0.5, 0, 0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.1) + ylab("blink response (arbitary V)") +
      scale_colour_discrete(labels=c("subject presses", "experimenter presses"), name="switching individual\npressing buttons")
    
    # faceted
    outside_facet <- ggplot(outside_df, aes(x=time, y=avg, group=blocks)) + geom_line(aes(colour=blocks), size=1) +
      theme_classic() + geom_vline(xintercept=0, linetype=2, alpha=0.4) + 
      scale_x_continuous(labels=function(x)x*1000, name = "time relative to button press (ms)", breaks=c(-0.5, 0, 0.5)) +
      geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se), alpha=0.3) + ylab("blink response (arbitary V)") +
      scale_colour_discrete(labels=c("subject presses", "experimenter presses"), name="switching individual\npressing buttons") + 
      facet_wrap(~blocks, ncol=2, strip.position = "left")
    
    if(grp=="control"){
      grouped_plot <- grouped_plot + scale_y_continuous(limits=c(min(blocks_df$avg - blocks_df$se)-0.01, max(blocks_df$avg + blocks_df$se)*5))
      facet_plot <- facet_plot + scale_y_continuous(limits=c(min(blocks_df$avg - blocks_df$se)-0.01, max(blocks_df$avg + blocks_df$se)*5))
      
      button_plot <- button_plot + scale_y_continuous(limits=c(min(button_df$avg - button_df$se)-0.01, max(button_df$avg + button_df$se)*5))
      
      switch_plot <- switch_plot + scale_y_continuous(limits=c(min(switch_df$avg - switch_df$se)-0.01, max(switch_df$avg + switch_df$se)*5))
      switch_facet <- switch_facet + scale_y_continuous(limits=c(min(switch_df$avg - switch_df$se)-0.01, max(switch_df$avg + switch_df$se)*5))
      
      outside_plot <- outside_plot + scale_y_continuous(limits=c(min(outside_df$avg - outside_df$se)-0.01, max(outside_df$avg + outside_df$se)*5))
      outside_facet <- outside_facet + scale_y_continuous(limits=c(min(outside_df$avg - outside_df$se)-0.01, max(outside_df$avg + outside_df$se)*5))
    }
    
    ggsave(grouped_plot, filename = "training - grouped plot.pdf", width=7, height=7*0.66)
    ggsave(facet_plot, filename = "training - individual blocks plot.pdf", width=7, height=7*0.66)
    ggsave(button_plot, filename = "buttons during training.pdf", width=7, height=7*0.66)
    ggsave(switch_plot, filename = "switching hands - grouped plot.pdf", width=7, height=7*0.66)
    ggsave(switch_facet, filename = "switching hands - separate plot.pdf", width=7, height=7*0.66)
    ggsave(outside_plot, filename = "outsider pressing - grouped.pdf", width=7, height=7*0.66)
    ggsave(outside_facet, filename = "outsider pressing - separate.pdf", width=7, height=7*0.66)
    
    setwd("..")
  }
  
  setwd("..")
}