#################################
# Figures for Report Appendix
# Project 3
# Rachel Weber
#################################

library(VIM)
library(ggplot2)

dat <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataRaw/vadata2.csv")
factor_vars <- c(2:4)
dat[,factor_vars] <- lapply(dat[,factor_vars] , factor)

dat <- dat[dat$proced != "2",]

table(dat$asa, dat$death30)
# there were very few people with asa of 1 and all survived beyond 30 days. Let's combine asa 1 and 2
dat[dat$asa == "1" & !is.na(dat$asa),]$asa <- "2"

# some BMIs are unrealistic. Let's remove those above 54
# there is a bmi of 2.5, which is crazy low, but realistically not impossible. We'll keep them for now
dat[dat$bmi > 60 & !is.na(dat$bmi),]$bmi <- NA

########### figure 1
pbox(dat[,c("bmi", "albumin")])


########### figure 2
est_data <- read.csv("C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_model.csv")

ggplot(est_data, aes(x = hospital, y = lower.ci)) + geom_point(color = "lightblue3") +
  geom_point(aes(x = hospital, y = upper.ci), color = "lightblue3") + 
  geom_point(aes(x = hospital, y = mean_pred), color = "coral3") +
  ggtitle("Estimated death rate for each hospital (period 39)") + theme_minimal() +
  ylab("Estimated Death Rate (%)") + xlab("Hospital")


no_alb <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_model_no_alb.csv")
bin_ci <- read.csv("C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/binom_ci.csv")


ggplot(bin_ci, aes(x = X)) + geom_segment(data = bin_ci,
                                                 mapping = aes(x = X, y = lower, xend = X, yend = upper), size = 1.5, color = "lightblue3") +
  geom_point(aes(x = X, y = Observed)) + theme_minimal() +
  geom_segment(data = est_data, mapping = aes(x = hospital + .3, y = lower.ci, xend = hospital + .3, yend = upper.ci), size = 1.5, color = "coral3", alpha = .5) +
  ggtitle("Observed and Expected Confidence Intervals for Hospital Rate in Period 39") + xlab("Hospital") + ylab("Death Rate (%)")



ggplot(bin_ci, aes(x = X)) + geom_segment(data = bin_ci,
                                          mapping = aes(x = X, y = lower, xend = X, yend = upper), size = 1.5, color = "lightblue3") +
  geom_point(aes(x = X, y = Observed)) + theme_minimal() +
  geom_segment(data = no_alb, mapping = aes(x = hospital + .3, y = lower.ci, xend = hospital + .3, yend = upper.ci), size = 1.5, color = "coral3", alpha = .5) +
  ggtitle("Observed and Expected Confidence Intervals for Hospital Rate in Period 39 (No Albumin)") + xlab("Hospital") + ylab("Death Rate (%)")
