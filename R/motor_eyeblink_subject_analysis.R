# subject dataframe plots
# comparisons of:
#      1. number of airpuffs by group
#      2. number of buttons
#      3. avg reaction time
#      4. avg ITI

# libraries
library(ggplot2)
library(dplyr)

# subject dataframe
setwd(choose.dir())
sub_df <- read.csv(file = "subject info.csv", header = T)

# boxplot of n airpuffs
airpuff_plot <- ggplot(sub_df, aes(x=group, y=n.airpuff, fill=group)) + geom_boxplot() + theme_classic() +
  scale_fill_manual(values = c("#F8766D", "#619CFF"), labels=c("pseudo-conditioning", "conditioning"), name="") +
  ylab("number of airpuffs during experiment") + scale_x_discrete(labels=c("pseudo-conditioning", "conditioning"))

ggsave(airpuff_plot, filename = "number of airpuff boxplot.pdf", width=7, height=7*0.66)

# stats for airpuff
sink("wilcoxon of number of airpuffs.txt")
wilcox.test(n.airpuff~group, data=sub_df)
sink()

# create a new dataframe for buttons
button_df <- as.data.frame(matrix(nrow=nrow(sub_df)*4, ncol=3))
colnames(button_df) <- c("group", "button", "n.presses")

# put in a data frame so each button for each subject has its own row
id_sub <- 1
for(i in seq(1, nrow(button_df), 4)){
  button_df$group[i] <- as.character(sub_df$group[id_sub])
  button_df$group[i + 1] <- as.character(sub_df$group[id_sub])
  button_df$group[i + 2] <- as.character(sub_df$group[id_sub])
  button_df$group[i + 3] <- as.character(sub_df$group[id_sub])
  
  button_df$button[i] <- "A"
  button_df$button[i + 1] <- "B"
  button_df$button[i + 2] <- "C"
  button_df$button[i + 3] <- "D"
  
  button_df$n.presses[i] <- sub_df$n.A[id_sub]
  button_df$n.presses[i + 1] <- sub_df$n.B[id_sub] + sub_df$n.BTest[id_sub]
  button_df$n.presses[i + 2] <- sub_df$n.C[id_sub] + sub_df$n.CTest[id_sub]
  button_df$n.presses[i + 3] <- sub_df$n.D[id_sub]
  
  id_sub <- id_sub + 1
}

# plot of the button
button_plot <- ggplot(button_df, aes(x=group, y=n.presses, fill=group)) + geom_boxplot() + facet_wrap(~button,ncol = 2, scales = "free_y") +
  theme_classic() + scale_fill_manual(values = c("#F8766D", "#619CFF"), labels=c("pseudo-conditioning", "conditioning"), name="") +
  ylab("number of button presses") + scale_x_discrete(labels=c("pseudo-conditioning", "conditioning")) + 
  theme(axis.text.x = element_text(angle = 355))

ggsave(button_plot, filename = "number of button presses.pdf", width=7, height=7*0.66)

# anova stats on button presses
button_anova <- aov(n.presses~group*button, data = button_df)

# anova stats in text file
sink("anova of button presses.txt")
summary(button_anova)
sink()

# average reaction time
RT_plot <- ggplot(sub_df, aes(x=group, y=avg.RT, fill=group)) + geom_boxplot() +
  theme_classic() + scale_fill_manual(values = c("#F8766D", "#619CFF"), labels=c("pseudo-conditioning", "conditioning"), name="") +
  ylab("average reaction time (ms)") + scale_x_discrete(labels=c("pseudo-conditioning", "conditioning")) + 
  scale_y_continuous(labels=function(x)x*1000)

# save RT plot
ggsave(RT_plot, filename = "average reaction time.pdf", width=7, height=7*0.66)

# rt plots stats
sink("reaction time stats.txt")
wilcox.test(avg.RT~group, data = sub_df)
sink()

# average ITI
ITI_plot <- RT_plot <- ggplot(sub_df, aes(x=group, y=avg.ITI, fill=group)) + geom_boxplot() +
  theme_classic() + scale_fill_manual(values = c("#F8766D", "#619CFF"), labels=c("pseudo-conditioning", "conditioning"), name="") +
  ylab("average inter-trial interval (s)") + scale_x_discrete(labels=c("pseudo-conditioning", "conditioning"))
ggsave(ITI_plot, filename = "intertrial interval.pdf", width = 7, height = 7*0.66)

# intertrial interval stats
sink("intertrial interval stats.txt")
wilcox.test(avg.ITI~group, data = sub_df)
sink()
