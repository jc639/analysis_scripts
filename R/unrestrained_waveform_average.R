# create average waveforms for all trial types for each subject, each day
library(dplyr)
library(ggplot2)

# Standard error function
se <- function(x) sqrt(var(x)/length(x))

# read in dataframe of waveforms
setwd(choose.dir())
trial_wave <- read.csv("trial waveforms.csv")

# subjects
subj <- unique(trial_wave$subject)

# group by subject, date, trial.type, time -> avg value and se
day_by_day <- trial_wave %>%
  group_by(subject, date, trial.type, time) %>%
  summarise(avg = mean(value),
            se = se(value))

# iterate through subjects and then dates
for(s in subj){
  
  # move into subjects folder
  setwd(as.character(s))
  subj_df <- filter(day_by_day, subject==s)
  max_value <- max(subj_df$avg + subj_df$se, na.rm = T) + 0.1
  min_value <- min(subj_df$avg - subj_df$se, na.rm = T) - 0.1
  
  dates <- unique(subj_df$date)
  
  dir.create("average waveforms")
  setwd("average waveforms")
  
  for(d in dates){
    date_df <- filter(subj_df, date == d)
    
    
    day_plot <- ggplot(date_df) + geom_line(aes(x=time, y=avg, colour=trial.type), size=1) + geom_ribbon(aes(x=time, ymax=avg+se, ymin=avg-se,
                                                                                                     group=trial.type),
                                                                                                 alpha=0.3) +
      theme_classic() + scale_x_continuous(limits = c(min(date_df$time, na.rm=T), max(date_df$time, na.rm = T)), breaks=c(-0.45, 0, 0.35),
                                           labels = function(x)x*1000, name="time relative to CS1 onset") +
      scale_y_continuous(limits = c(min_value, max_value), name = "average value") +
      scale_color_manual(values = c("#C77CFF", "#F8766D", "#00BA38","#619CFF"), guide=guide_legend(title = "trial type"))
    
    ggsave(day_plot, filename = paste0(s, " ", d, ".pdf"), width = 9, height = 6)
  }
  setwd("..")
  setwd("..")
}
