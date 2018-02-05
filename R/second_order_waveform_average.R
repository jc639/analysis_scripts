# daily average waveforms for each subject by session
# need to correct session from other csv to block congruently
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

# standard error function
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
} 

# choose folder with each subject folder in
overall_dir <- choose.dir()
setwd(overall_dir)
folder_list <- list.dirs(recursive = F)

# iterate through the subject folders
for (d in folder_list){
  setwd(d)
  
  # create folder for animation
  dir.create("animation")
  
  # list files
  file_list <- list.files()
  
  # read in waveform csv, but also need response to write correct session
  trial_waveform <- read.csv(file_list[grep(x = file_list, pattern = "waveform")], header = T)
  trial_response <- read.csv(file_list[grep(x = file_list, pattern = "response")], header = T)
  
  if(d!="./Ben"){
    # generate session column for trial waveform, first concatenate date and session n
    trial_waveform$session <- paste0(trial_waveform$session.n, trial_waveform$date)
    
    # fill in the correct session, find unique combination of session.n, date and session in response df
    session_response_df <- unique(cbind(trial_response$session.n, trial_response$date, trial_response$session))
    # concatenate these unique identifiers
    session_response_df[, 1] <- paste0(session_response_df[,1], session_response_df[,2])
    # make column 2 the session number
    session_response_df[, 2] <- session_response_df[,3]
    # drop the third column
    session_response_df <- session_response_df[,1:2]
    
    # vectorised which function, identify which session_response[,1] is equivalent to the given trial_waveform$session
    # and then uses that to return the numeric session number
    trial_waveform$session <- as.numeric(session_response_df[,2][unname(sapply(trial_waveform$session,
                                                                               FUN = function(x) which(session_response_df[,1]==x)))])
  } else {
    uni_date <- unique(trial_waveform$date)
    trial_waveform$session <- unname(sapply(trial_waveform$date, FUN = function(x) which(uni_date==x)))
    trial_waveform <- filter(trial_waveform, trial.type!="CS1")
  }
  
  # filter so we are only working on CR/No CR trials
  trial_waveform <- filter(trial_waveform, response%in%c("CR", "NO CR"))
  
  # iterate through and create plots for each day
  # first create time line diagram to grid arrange above the plot
  onset_diagram <- ggplot()+geom_rect(aes(xmin = -0.35, xmax=0.06, ymin=-1.2, ymax=-0.2), fill="#619CFF")+
    geom_rect(aes(xmin=0, xmax=0.42, ymax=-1.2, ymin=-2.2), fill= "#F8766D") +
    geom_rect(aes(xmin=0.35, xmax=0.42, ymin=-3.3, ymax=-2.2), fill="grey") +
    theme_bw() + theme(plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border= element_blank(),
                       axis.line.x = element_line(color="black"),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())+
    geom_text(aes(x=-0.14, y=-0.7,  label="CS2"), size=6) +
    geom_text(aes(x=0.21, y=-1.7, label="CS1"), size=6) +
    geom_text(aes(x=0.385, y=-2.7, label="US"), size=6) +
    xlab("time relative to CS1 stimulus onset (ms)") +
    scale_x_continuous(limits = c(-0.5, 0.45), breaks = c(-0.35, 0, 0.35), labels = function(x)x*1000, position = "top",
                       expand=c(0,0))+scale_y_continuous(expand=c(0,0))
    
  
  # vector of unique sessions
  uni_session <- unique(trial_waveform$session)
  
  # filter through the sessions, creating a plot of average responses for each trial type
  sess_df <- trial_waveform %>%
      group_by(session, trial.type, time) %>%
      summarise(avg = mean(value),
                se = se(value))
  
  # iterate through the sessions and filter generate folder, save two plots
  for(sess in uni_session){
    sess_plot <- ggplot(filter(sess_df, session==sess)) + geom_line(aes(x=time, y=avg, colour=trial.type), size=1) +
      geom_ribbon(aes(x=time, ymin=avg-se, ymax=avg+se, group=trial.type), alpha=0.3) + 
      theme_classic() + scale_x_continuous(limits = c(min(sess_df$time), max(sess_df$time)), breaks = c(-350, 0, 350), expand = c(0,0)) +
      xlab("time relative to CS1 stimulus onset (ms)") + ylab("blink size (mm)") +
      scale_colour_manual(values = c("#F8766D", "#00BA38","#619CFF"), guide=guide_legend(title = "trial type:")) +
      scale_y_continuous(limits = c(min(sess_df$avg)-.2, max(sess_df$avg)+.2)) + theme(legend.position = "top")
    
    grid_plot <- arrangeGrob(onset_diagram, sess_plot, ncol=1, heights=unit(c(1.5, 6), c("in", "in")), widths=unit(9, "in"))
    
    # create folder for each session
    dir.create(paste0("session ", sess, " plots"))
    setwd(paste0("session ", sess, " plots"))
    
    # save plots
    ggsave(grid_plot, filename = paste0("session ", sess, " plot with onset diagram.pdf"), width=9,
           height=7.5)
    
    # add two the plot
    sess_plot2 <- sess_plot + geom_vline(xintercept = c(-350, 0, 350), alpha=0.4, linetype=2) +
      annotate("text", x = -290, y = max(sess_df$avg)+.2, label="CS2 onset") +
      annotate("text", x = 60, y = max(sess_df$avg)+.2, label="CS1 onset") +
      annotate("text", x = 400, y = max(sess_df$avg)+.2, label="US onset")
    
    ggsave(sess_plot2, filename = paste0("session ", sess, " plot.pdf"), width=9, height=6)
    
    # move back to subject folder
    setwd("..")
    
    # move into animation folder - this is to generate a gif of plots
    setwd("animation")
    
    # add a title so can be identified in animation
    ani_plot <- sess_plot + ggtitle(paste0("session ", sess))
    ggsave(ani_plot, filename = paste0("session ", sess, ".pdf"), width=6, height=4)
    setwd("..")
  }  
  setwd("..")
}

