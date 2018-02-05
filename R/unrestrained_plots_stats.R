# unrestrained analysis
# so two groups:
#     light primary:
#       - CS1US - flashing light 350ms CS-US interval, US = airpuff
#       - CS2CS1 - 450ms interval 1kHz tone to CS1 light - 800ms in total
#
#     sound primary:
#       - CS1US - sound 1kHz tone 350ms interval
#       - CS2CS1 - 450ms interval Light
#
# Output of this script:
#   - Individual learning curves for each subject
#   - onsets by session
#   - group by primary CS
#   - stats

# libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)

# standard error function
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
} 

# select file location
folder_dir <- setwd(choose.dir())

# read in csv of trials
trial_response <- read.csv("response dataframe.csv", header = T)

# add group to the dataframe (light/sound)
sound_group <- c("Bashful", "Doc", "Sleepy")
trial_response$group <- ifelse(trial_response$subject %in% sound_group, "sound", "light")

# create session numbers from dates
uni_date <- unique(trial_response$date)
trial_response$session <- sapply(trial_response$date, FUN = function(x) which(uni_date==x))

# create individual learning graphs
cr_df <- trial_response %>%
  filter(trial.type!="CS1", response %in% c("CR", "NO CR")) %>%
  group_by(subject, session, trial.type) %>%
  summarise(perc=(length(response[response=="CR"])/n()) * 100,
            group = unique(group),
            first.order = unique(first.order))

# block the sessions to make it bit more
cr_df$cut <- cut(cr_df$session, breaks = c(0, 3, 6, 9, 12, 15, 18, 20, 22, 25, 28, 31, 34, 37, 40))
unique_cut <- unique(cr_df$cut)
cr_df$block <- sapply(cr_df$cut, FUN = function(x) which(unique_cut == x))


# need to add NAs for CS2CS1 and CS1US for blocks where CS1US was only given
uni_sub <- as.character(unique(cr_df$subject))
cs1_only <- unique(cr_df$block[cr_df$first.order==T])
tmp_df <- data.frame(subject = rep(uni_sub, 2), session = rep(NA, 12),
                     trial.type = c(rep("CS2CS1", 6), rep("CS2", 6)), perc = rep(NA, 12),
                     group = ifelse(rep(uni_sub, 2) %in% sound_group, "sound", "light"),
                     first.order = rep(T, 12), cut = rep(NA, 12), block = c(rep(cs1_only[1], 6), rep(cs1_only[2], 6)))

# need to convert to data frame
cr_df <- rbind(as.data.frame(cr_df), tmp_df)

# df for plots of cr incidence
cr_plot_df <- cr_df %>%
  group_by(subject, block, trial.type) %>%
  summarise(avg = mean(perc),
            se = se(perc))

# create blocks in trial response so onset plots can be done on the same
trial_response$cut <- cut(trial_response$session, breaks = c(0, 3, 6, 9, 12, 15, 18, 20, 22, 25, 28, 31, 34, 37, 40))
unique_cut <- unique(cr_df$cut)
trial_response$block <- sapply(trial_response$cut, FUN = function(x) which(unique_cut == x))

# create individual onset graphs
onset_df <- trial_response %>%
  filter(trial.type!="CS1", response %in% c("CR"), !is.infinite(onset), trial.type%in%c("CS1US", "CS2CS1")) %>%
  group_by(subject, block, trial.type) %>%
  summarise(avg_onset = mean(onset, na.rm=T),
            se = se(onset))

# use temp df change names to bind with onset df
tmp_df <- tmp_df[, c("subject", "block", "perc", "cut", "trial.type")]
colnames(tmp_df)[which(colnames(tmp_df)=="perc")] <- "avg_onset" 
colnames(tmp_df)[which(colnames(tmp_df)=="cut")] <- "se"

# onset df
onset_df <- rbind(as.data.frame(onset_df), tmp_df)
# CS2 is in tmp_df, don't need for onsets as so few CRS
onset_df <- filter(onset_df, trial.type!="CS2")

# unique subjects, create a folder for each
uni_sub <- as.character(unique(trial_response$subject))
sapply(uni_sub, FUN = function(x) dir.create(x))

# create plots, first stats/modelling next
for(s in uni_sub){
  setwd(s)
  
  # filter trial response to just subject
  subj_response <- filter(trial_response, subject==s, trial.type%in%c("CS1US", "CS2CS1"),
                          response=="CR")
  
  # create blocked CR plot for subject
  cr_plot <- ggplot(filter(cr_plot_df, subject==s), aes(x=block, y=avg, colour=trial.type)) + geom_point() + geom_line() +
    geom_errorbar(aes(x=block, ymin=avg-se, ymax=avg+se), width=0.25) + theme_classic() + scale_y_continuous(limits=c(0, 100)) +
    ylab("percentage of CRs") + scale_colour_manual(values = c("#F8766D", "#00BA38","#619CFF"), guide=guide_legend(title = "trial type")) + 
    annotate("rect", xmin=6.5, xmax=8.5, ymin=-Inf, ymax=Inf, alpha=0.2) + annotate("text", x=7.5,y=95, label="first order\n\tonly")
  
  
  # cs2 probabilty of CR by session, for predicted plot
  cs2_prob <- trial_response %>%
    filter(subject == s, trial.type=="CS2", response%in%c("CR", "NO CR")) %>%
    group_by(session) %>%
    summarise(total=n(),
              CRs = length(response[response=="CR"]),
              prob = (CRs/total)*100)
  
  # create a df for the glm, need to indicate CR with 1
  cs2_cr <- trial_response %>%
    filter(subject == s, trial.type=="CS2", response%in%c("CR", "NO CR"))
  cs2_cr$CR <- ifelse(cs2_cr$response=="CR", 1, 0)
  
  # need to lower session numbers to account for the first order period
  cs2_cr$session[cs2_cr$session >= 23] <- cs2_cr$session[cs2_cr$session >= 23] - 4
  
  # create session glm
  cs2_sess_glm <- glm(CR~session, data=cs2_cr, family = binomial(link = "logit"))
  
  # create a prediction df
  pred_df <- data.frame(session = 1:max(cs2_cr$session), 
                        pred = predict(cs2_sess_glm, newdata = data.frame(session=1:max(cs2_cr$session)),
                                       type = "response"))
  pred_df$pred <- pred_df$pred * 100
  
  # add back the missing sessions, and add NAs for the gap
  pred_df$session[pred_df$session >= 19] <- pred_df$session[pred_df$session >= 19] + 4
  gap_df <- data.frame(session=c(19,20,21,22), pred=rep(NA, 4))
  pred_df <- rbind(pred_df, gap_df)
  
  # session probabilty plot
  cr_sess_prob <- ggplot() + geom_point(data=cs2_prob, aes(x=session, y=prob), colour="#00BA38", size=1.5) +
    geom_line(data=pred_df, aes(x=session, y=pred), colour="#00BA38", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(pred_df$session)), breaks=seq(5, max(pred_df$session), 5)) +
    ylab("Probability of CR") + xlab("session number") + scale_y_continuous(limits = c(0,100)) +
    annotate("rect", xmin=18.5, xmax=22.5, ymin=-Inf, ymax=Inf, alpha=0.2) + annotate("text", x=20.5,y=50, label="first\norder\nonly")
  
  
  
  ##################### same but for trial ###################################
  # cs2 probability for trial
  cs2_prob <- trial_response %>%
    filter(subject == s, trial.type=="CS2", response%in%c("CR", "NO CR")) %>%
    group_by(trial.n) %>%
    summarise(total=n(),
              CRs = length(response[response=="CR"]),
              prob = (CRs/total)*100)
  
  # cs2 generalised linear model for trial.n predictor
  cs2_trial_glm <- glm(CR~trial.n, data=cs2_cr, family = binomial(link = "logit"))
  
  # prediction dataframe
  pred_df <- data.frame(trial.n = 1:max(cs2_cr$trial.n), 
                        pred = predict(cs2_trial_glm, newdata = data.frame(trial.n=1:max(cs2_cr$trial.n)),
                                       type = "response"))
  pred_df$pred <- pred_df$pred * 100
  
  # plot of trial predicted against observed 
  cr_trial_prob <- ggplot() + geom_point(data=cs2_prob, aes(x=trial.n, y=prob), colour="#00BA38", size=1.5) +
    geom_line(data=pred_df, aes(x=trial.n, y=pred), colour="#00BA38", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(pred_df$trial.n)), breaks=seq(1, max(pred_df$trial.n), 1)) +
    ylab("Probability of CR") + xlab("trial number") + scale_y_continuous(limits = c(0,100))
  
  # create both session and trial.n glm and interaction
  simple_cs2_glm <- glm(CR~session+trial.n, data=cs2_cr, family = binomial(link = "logit"))
  int_cs2_glm <- glm(CR~session*trial.n, data=cs2_cr, family = binomial(link = "logit"))
  
  # sink the CS2 GLMs into a text file
  sink(paste0(s, " - CS2 response glms.txt"))
  cat("SESSION GLM\n")
  print(summary(cs2_sess_glm))
  cat("\n\n\n\n\n")
  cat("TRIAL GLM\n")
  print(summary(cs2_trial_glm))
  cat("\n\n\n\n\n")
  cat("BOTH VARIABLE GLM\n")
  print(summary(simple_cs2_glm))
  cat("\n\n\n\n\n")
  cat("INTERACTION BETWEEN VARIABLES\n")
  print(summary(int_cs2_glm))
  sink()
  
  # onset blocked by the same parameters as CR plot
  onset_plot <- ggplot(filter(onset_df, subject==s), aes(x=block, y=avg_onset, colour=trial.type)) + geom_point() + geom_line() +
    geom_errorbar(aes(x=block, ymin=avg_onset-se, ymax=avg_onset+se), width=0.25) + theme_classic() + 
    scale_y_continuous(limits = c(-0.45, 0.35), breaks = c(-0.45, 0, 0.35), labels = function(x)x*1000) +
    geom_hline(aes(yintercept=0), linetype=2, alpha=0.4) + annotate("rect", xmin=6.5, xmax=8.5, ymin=-Inf, ymax=Inf, alpha=0.2) + 
    annotate("text", x=7.5,y=0.3, label="first order\n\tonly") + ylab("average onset (ms)") +
    scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type"))
  
  # onset dot plot for sessions, CS1US
  cs1_session_plot <- ggplot(filter(subj_response, trial.type=="CS1US", onset < 0.35), aes(x=session, y=onset)) + geom_point(colour="#F8766D", alpha=0.5) + 
    theme_classic() + scale_y_continuous(limits = c(-0.46,0.36), labels = function(y)y*1000, breaks = c(-0.45, 0, 0.35)) +
    scale_x_continuous(limits = c(1, max(subj_response$session)), breaks=seq(5, max(subj_response$session), 5)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS1-US")
  
  # onset dot plot for session, CS2CS1
  cs2_session_plot <- ggplot(filter(subj_response, trial.type=="CS2CS1", onset < 0.35), aes(x=session, y=onset)) + geom_point(colour="#619CFF", alpha=0.5) + 
    theme_classic() + scale_y_continuous(limits = c(-0.46,0.36), labels = function(y)y*1000, breaks = c(-0.45, 0, 0.35)) +
    scale_x_continuous(limits = c(1, max(subj_response$session)), breaks=seq(5, max(subj_response$session), 5)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS2-CS1") + annotate("rect", xmin=18.5, xmax=22.5, ymin=-Inf, ymax=Inf, alpha=0.2) + 
    annotate("text", x=20.5,y=0.2, label="first\norder\nonly")
  
  # early indicator - 1 if onset < 0 
  subj_response$early <- ifelse(subj_response$onset < 0, 1, 0)
  
  # probability of early response for CS2CS1 trials
  prob_df <- subj_response %>%
    filter(trial.type=="CS2CS1") %>%
    group_by(session) %>%
    summarise(total=n(),
              early=length(onset[onset < 0]),
              prob = early/total)
  
  
  # create an early onset glm from df, but need to account for no sessions during first order (-4 after session 23)
  early_onset <- subj_response %>%
    filter(trial.type=="CS2CS1")
  early_onset$session[early_onset$session >= 23] <- early_onset$session[early_onset$session >= 23] - 4
  
  # glm
  early_sess_onset_glm <- glm(early~session, data = early_onset, family = binomial(link = "logit"))
  
  # create a prediction df
  pred_df <- data.frame(session = 1:max(early_onset$session), 
                        pred = predict(early_sess_onset_glm, newdata = data.frame(session=1:max(early_onset$session) + 4),
                                                                             type = "response"))
  # add back the missing sessions, and add NAs for the gap
  pred_df$session[pred_df$session >= 19] <- pred_df$session[pred_df$session >= 19] + 4
  gap_df <- data.frame(session=c(19,20,21,22), pred=rep(NA, 4))
  pred_df <- rbind(pred_df, gap_df)
  
  # session probabilty plot
  sess_prob_plot <- ggplot() + geom_point(data=prob_df, aes(x=session, y=prob), colour="#619CFF", size=1.5) +
    geom_line(data=pred_df, aes(x=session, y=pred), colour="#619CFF", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(pred_df$session)), breaks=seq(5, max(pred_df$session), 5)) +
    ylab("Probability of early onset (< 0 ms)") + xlab("session number") + scale_y_continuous(limits = c(0,1)) +
    annotate("rect", xmin=18.5, xmax=22.5, ymin=-Inf, ymax=Inf, alpha=0.2) + annotate("text", x=20.5,y=0.5, label="first\norder\nonly")
  
  
  ##################### same but for trial ###################################
  
  
  # onset dot plot for trials, CS1US
  cs1_trial_plot <- ggplot(filter(subj_response, trial.type=="CS1US", trial.n<=30, onset < 0.35), aes(x=trial.n, y=onset)) + geom_point(colour="#F8766D", alpha=0.5) + 
    theme_classic() + scale_y_continuous(limits = c(-0.46,0.36), labels = function(y)y*1000, breaks = c(-0.45, 0, 0.35)) +
    scale_x_continuous(limits = c(1, 30), breaks=seq(5, 30, 5)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS1-US")
  
  # onset dot plot for session, CS2CS1
  cs2_trial_plot <- ggplot(filter(subj_response, trial.type=="CS2CS1", onset < 0.35), aes(x=trial.n, y=onset)) + geom_point(colour="#619CFF", alpha=0.5) + 
    theme_classic() + scale_y_continuous(limits = c(-0.46,0.36), labels = function(y)y*1000, breaks = c(-0.45, 0, 0.35)) +
    scale_x_continuous(limits = c(1, 30), breaks=seq(5, 30, 5)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS2-CS1")
  
  
  # probability of early response for CS2CS1 trials
  prob_df <- subj_response %>%
    filter(trial.type=="CS2CS1") %>%
    group_by(trial.n) %>%
    summarise(total=n(),
              early=length(onset[onset < 0]),
              prob = early/total)
  
  
  # create an early onset glm from df, but need to account for no sessions during first order (-4 after session 23)
  # glm
  early_trial_onset_glm <- glm(early~trial.n, data = early_onset, family = binomial(link = "logit"))
  
  # create a prediction df
  pred_df <- data.frame(trial.n = 1:max(early_onset$trial.n), 
                        pred = predict(early_trial_onset_glm, newdata = data.frame(trial.n=1:max(early_onset$trial.n)),
                                       type = "response"))
 
  # session probabilty plot
  trial_prob_plot <- ggplot() + geom_point(data=prob_df, aes(x=trial.n, y=prob), colour="#619CFF", size=1.5) +
    geom_line(data=pred_df, aes(x=trial.n, y=pred), colour="#619CFF", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(early_onset$trial.n)), breaks=seq(5, max(early_onset$trial.n), 5)) +
    ylab("Probability of early onset (< 0 ms)") + xlab("trial number with session") + scale_y_continuous(limits = c(0,1))
  
  
  # glm of both variables and interaction
  early_simple <- glm(early~session+trial.n, data=early_onset, family = binomial(link = "logit"))
  early_int <- glm(early~session*trial.n, data=early_onset, family = binomial(link = "logit"))
  
  
  # sink the early onset glms in to a text file
  sink(paste0(s, " - glm on early onset in CS2CS1.txt"))
  cat("SESSION GLM\n")
  print(summary(early_sess_onset_glm))
  cat("\n\n\n\n\n")
  cat("TRIAL GLM\n")
  print(summary(early_trial_onset_glm))
  cat("\n\n\n\n\n")
  cat("BOTH VARIABLES GLM\n")
  print(summary(early_simple))
  cat("\n\n\n\n\n")
  print(summary(early_int))
  sink()
  
  # make plots of CS1US and CS2CS1 probabilities against predicted
  primary_probs <- trial_response %>%
    filter(subject==s, trial.type %in% c("CS1US", "CS2CS1"), response %in% c("CR", "NO CR")) %>%
    group_by(session, trial.type) %>%
    summarise(prob = (length(response[response=="CR"])/n())*100)
  
  # make a sub df to pass to glm 
  sub_df <- trial_response %>%
    filter(subject == s, trial.type%in%c("CS1US", "CS2CS1"), response%in%c("CR", "NO CR"))
  sub_df$success <- ifelse(sub_df$response=="CR", 1, 0)
  
  #glms
  simple_primary_glm <- glm(success~session+trial.type, data=sub_df, family = binomial(link = "logit"))
  int_primary_glm <- glm(success~session*trial.type, data=sub_df, family = binomial(link = "logit"))
  
  # sink the results of these glms into a text file
  sink(paste0(s, " - glm of primary learning.txt"))
  cat("SIMPLE GLM\n")
  print(summary(simple_primary_glm))
  cat("\n\n\n\n\n\n")
  cat("INTERACTION GLM\n")
  print(summary(int_primary_glm))
  sink()
  
  # predictions
  primary_probs$simple <- predict(simple_primary_glm, newdata = primary_probs[, 1:2], type = "response") * 100
  primary_probs$int <- predict(int_primary_glm, newdata = primary_probs[, 1:2], type = "response") * 100
  
  # need to add on NAs for CS2C1
  gap_df <- data.frame(session=c(19,20,21,22), trial.type=rep("CS2CS1", 4), prob = rep(NA, 4),
                       simple = rep(NA, 4), int = rep(NA, 4))
  primary_probs <- rbind(as.data.frame(primary_probs), gap_df)
  
  # make two plots, one for simple, one for interaction glm
  simple_prim <- ggplot(primary_probs) + geom_point(aes(x=session, y=prob, colour=trial.type), alpha=0.6) + 
    geom_line(aes(x=session, y=simple, colour=trial.type), size=1) + theme_classic() +
    scale_y_continuous(limits = c(0,100), name = "percentage of CRs") + 
    scale_color_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type")) +
    scale_x_continuous(limits=c(0.5,max(primary_probs$session)+0.5), breaks=seq(10, max(primary_probs$session), 10))
  
  int_prim <- ggplot(primary_probs) + geom_point(aes(x=session, y=prob, colour=trial.type), alpha=0.6) + 
    geom_line(aes(x=session, y=int, colour=trial.type), size=1) + theme_classic() + 
    scale_y_continuous(limits = c(0,100), name = "percentage of CRs") + 
    scale_color_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type")) +
    scale_x_continuous(limits=c(0.5,max(primary_probs$session)+0.5), breaks=seq(10, max(primary_probs$session), 10))
  
  # make a sub df to pass to chisq 
  sub_df <- trial_response %>%
    filter(subject == s, trial.type%in%c("CS1US", "CS2CS1"))
  
  # Chisq of early onset, is there an overall difference in CS1US/CS2CS1 early onset
  cs1_early <- ifelse(sub_df$trial.type=="CS1US" & sub_df$response=="not flat baseline", 1, 0)
  cs2_early <- ifelse(sub_df$trial.type=="CS2CS1" & !is.na(sub_df$onset) &
                        sub_df$onset < 0, 1, 0)
  cs1_tots <- length(cs1_early) - sum(cs1_early)
  cs2_tots <- length(cs2_early) - sum(cs2_early)
  
  chi_mat <- as.table(rbind(c(cs1_tots, cs2_tots), c(sum(cs1_early), sum(cs2_early))))
  dimnames(chi_mat) <- list(timing=c("Not Early", "Early"),
                          trial=c("CS1-US", "CS2-CS1"))
  
  # chisq matrix
  chi_sq <- chisq.test(chi_mat)
  chi_table <- tableGrob(chi_mat)
  
  # save table
  ggsave(chi_table, filename = paste0(s, " - chi table.pdf"))
  
  sink(paste0(s, " - results of chi sq.txt"))
  print(chi_sq)
  sink()
  
  
  # save all plots
  grid_plot <- arrangeGrob(cs1_session_plot, cs2_session_plot, sess_prob_plot, 
                           ncol=1, heights=unit(c(3,3,3), rep("in", 4)), widths=unit(6, "in"))
  
  ggsave(grid_plot, filename = paste0(s, " - session onsets.pdf")  , height = 12, width = 8, units = "in")
  
  grid_plot <- arrangeGrob(cs1_trial_plot, cs2_trial_plot, trial_prob_plot, 
                           ncol=1, heights=unit(c(3,3,3), rep("in", 4)), widths=unit(6, "in"))
  ggsave(grid_plot, filename = paste0(s, " - trial onsets.pdf")  , height = 12, width = 8, units = "in")
  
  ggsave(cr_plot, filename = paste0(s, " - cr percentages.pdf"), width=9, height=6)
  
  ggsave(cr_sess_prob, filename = paste0(s, " - probability of CS2 by session.pdf"), width = 6, height = 4)
  
  ggsave(cr_trial_prob, filename = paste0(s, " - probability of CS2 by trial.pdf"), width = 6, height = 4)
  
  ggsave(simple_prim, filename = paste0(s, " - predicted probabilities for primary learning - simple model.pdf"),
         width = 6, height = 4)
  
  ggsave(int_prim, filename = paste0(s, " - predicted probabilities for primary learning - interaction model.pdf"),
         width = 6, height = 4)
  
  ggsave(onset_plot, filename = paste0(s, " - onset blocked.pdf"), width=9, height=6)
  
  
  setwd("..")
}


