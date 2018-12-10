##############################
# Exploratory
# Rachel Weber
# Project 3: Missing Data
#############################

library(ggplot2)
library(VIM)
library(gridExtra)


dat <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataRaw/vadata2.csv")

factor_vars <- c(1:3,9)
dat[,factor_vars] <- lapply(dat[,factor_vars] , factor)

sapply(dat, function(x) sum(is.na(x)))
  # lots of missing albumin 


aggr(dat, numbers = T, prop = c(T, F))
  # well clearly we've got an albumin problem

dat$miss_a <- ifelse(is.na(dat$albumin), 1, 0)
dat$miss_a <- factor(dat$miss_a)



ggplot(dat, aes(x = weight, y = albumin, color = sixmonth)) + geom_point()

ggplot(dat, aes(x = height, y = albumin)) + geom_point()

ggplot(dat, aes(x = bmi, y = albumin)) + geom_point()
summary(dat$bmi)
miss <- dat[is.na(dat$albumin),]
summary(miss$bmi)

ggplot(dat, aes(x = bmi, fill = miss_a)) +
  geom_histogram(alpha = .5, position = 'identity')
  # perfect overlap...


# maybe albumin is missing in certain categories...
spineMiss(dat[,c("hospcode", "albumin")])
spineMiss(dat[,c("sixmonth", "albumin")])
spineMiss(dat[,c("death30", "albumin")])
spineMiss(dat[,c("proced", "albumin")])
spineMiss(dat[,c("asa", "albumin")])


# boxplot for missing vs. not
pbox(dat[,c("weight", "albumin")])
pbox(dat[,c("height", "albumin")])
pbox(dat[,c("bmi", "albumin")])

nrow(dat[dat$death30 == "1",])/nrow(dat) # .04
nrow(miss[miss$death30 == "1",])/nrow(miss) # .06
  # nearly identical

# I think it's missing at random...............