##################################
# Project 0 - BIOS 6624
# Rachel Weber
# Exploratory Analysis
# Edited: 8/28/18
##################################

library(ggplot2)
library(gridExtra)
library(car)
library(dplyr)

# read data
dat <- read.csv(file = "C:/Users/weberra/Downloads/Project0_dental_data.csv", header = T)

# make treatment factors
dat$trtgroup <- as.factor(dat$trtgroup)

# round variables to 2 sig figs
dat[,c(8:11)] <- round(dat[,c(8:11)], 2)

# attachment by treatment
g1 <- ggplot(dat, aes(x = trtgroup, y = attachbase)) + geom_point()
g2 <- ggplot(dat, aes(x = trtgroup, y = attach1year)) + geom_point()

grid.arrange(g1, g2)

# pd by treatment
g3 <- ggplot(dat, aes(x = trtgroup, y = pdbase)) + geom_point() + ylim(2,5) +
g4 <- ggplot(dat, aes(x = trtgroup, y = pd1year)) + geom_point() + ylim(2,5) +
  
grid.arrange(g3, g4)

# calculate changes
dat$delta_attach <- dat$attachbase - dat$attach1year
dat$delta_depth <- dat$pdbase - dat$pd1year

# plot changes
ggplot(dat, aes(attachbase, attach1year)) + geom_point()
ggplot(dat, aes(pdbase, pd1year)) + geom_point()
ggplot(dat, aes(trtgroup, delta_attach)) + geom_point()
ggplot(dat, aes(trtgroup, delta_depth)) + geom_point()

# mean attachment change by group
dat %>% 
  group_by(trtgroup) %>% 
  summarize(mean_delta_attach = mean(delta_attach, na.rm = T))

# mean depth change by group
dat %>% 
  group_by(trtgroup) %>% 
  summarize(mean_delta_depth = mean(delta_depth, na.rm = T))

################# linear models ####################
# get complete cases so we can do predictions and R-squared
dat1 <- dat[complete.cases(dat),]

summary(lm(delta_attach ~ trtgroup, data = dat))
plot(lm(delta_attach ~ trtgroup, data = dat))
  # QQ plot shows underprediction at low and high...may not be linear relationship
# model 1: change in attachment adjusted by baseline
m1 <- glm(delta_attach ~ trtgroup + attachbase, data = dat1)
Anova(m1)
cor(predict(m1), dat1$delta_attach)
  #.473
m3 <- glm(delta_attach ~ trtgroup + trtgroup^2 + attachbase, data = dat1)
Anova(m3)
cor(predict(m3), dat1$delta_attach)
#.473

# if we adjust for other demographic variables...
m2 <- glm(delta_attach ~ trtgroup + attachbase + smoker + sex + race + age, data = dat1)
  #smoking is a nonsignificant predictor
Anova(m2)
cor(predict(m2), dat1$delta_attach)
  #.525


