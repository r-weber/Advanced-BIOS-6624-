##################################
# Exploratory Analysis
# Author: Rachel Weber
# BIOS 6624 Project 1
# Created: Sept 14, 2018
# Edited:
##################################

# QOI: We are interested in understanding how treatment response 2 years after 
# initiating HAART differs between subjects who report using hard drugs, such as heroine and 
# cocaine, at baseline and other subjects, who did not report hard drug use at baseline.

library(CIDAtools)
library(ggplot2)
library(data.table)
library(car)
library(magrittr)
library(tidyverse)
library(rjags)
library(mcmcse)

dat <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataRaw/hiv_6624_final.csv", sep = ",")

setnames(dat, tolower(names(dat)))

bin_col <- c(5:7, 9:14, 18, 20:23, 26:29, 31:32, 34)
dat[, bin_col] <- lapply(dat[, bin_col], as.factor)

# table 1

# get first record
first <- dat[dat$years == 0,]
Table1(first, c("race", "smoke", "age", "bmi", "adh", "income", "hard_drugs"), colvar = NULL, verbose = T, incl_missing = F)

# missing values
sapply(dat, function(x) sum(is.na(x)))
  # triglycerides, LDL, tchol all above 900 missing
  # leu3n, vload, adh, bmi above 200 missing

############################ plotting outcomes ###############################
# plotting vload over time
# subset down to 2 years and less, since that's the QOI
dat1 <- dat[dat$years <= 2,]
# dat1 <- dat1[!is.na(dat1$hard_drugs),]

ggplot(dat1, aes(years, vload, group = newid)) + 
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) +
  facet_grid(.~ hard_drugs) + ylim(0, 20000)
  # max vload is 190 mil... shrunk graph to show comparable dist'n of means
  # both groups show improvement with time-esp from 0 to 1 year

# T cells over time
ggplot(dat1, aes(years, leu3n, group = newid)) + 
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) +
  facet_grid(.~ hard_drugs)
  # both show nearly identical improvement over time

# physical QOL over time
ggplot(dat1, aes(years, agg_phys, group = newid)) + 
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) +
  facet_grid(.~ hard_drugs)
  # trend differently; people on drugs report improvement over time
  # those not using hard drugs have decline from year 0

# mental QOL over time
ggplot(dat1, aes(years, agg_ment, group = newid)) + 
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) +
  facet_grid(.~ hard_drugs)
  # odd trend for people on drugs; ultimately improvement at 2 years though
  # never reach same mental health level as non-users

####################### data management #############################

# transform by log base 10
dat1$log_vload <- log10(dat1$vload)
dat1$log_tcell <- log10(dat1$leu3n)

# noticed 2 people with BMIs of 999; making those values--and all clearly unrealistic--NA
dat1[!is.na(dat1$bmi) & dat1$bmi > 54,]$bmi <- NA
dat1[!is.na(dat1$bmi) & dat1$bmi < 0,]$bmi <- NA

write.table(dat1, "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataProcessed/hiv_for_graphs.csv", sep = ",")

# reshape long to wide
dat2 <- reshape(dat1, idvar = "newid", timevar = "years", direction = "wide")

# make change outcomes
dat2$delta_vload <-  dat2$log_vload.0-dat2$log_vload.2
dat2$delta_tcell <- dat2$log_tcell.2-dat2$log_tcell.0
dat2$delta_mental <- dat2$agg_ment.2-dat2$agg_ment.0
dat2$delta_phys <- dat2$agg_phys.2-dat2$agg_phys.0

# drop columns from year 1
dat2 <- dat2[, -grep("1$", colnames(dat2))]
dat2 <- dat2[, -grep("^x", colnames(dat2))]

# condense adherence categories
levels(dat2$adh.2) <- c("0", "0", "1", "1")

# condense smoking 
levels(dat2$smoke.0) <- c("0", "0", "1")
levels(dat2$smoke.2) <- c("0", "0", "1")

# condense race
levels(dat2$race.0) <- c("White", "White", "Other", "Other", "Other", "Other")
levels(dat2$race.2) <- c("White", "White", "Other", "Other", "Other", "Other")

# condense education
levels(dat2$educbas.0) <- c("< College", "< College", "< College", "< College", 
                            "4 years/Degree", "> College", "> College")

write.table(dat2, "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataProcessed/hiv_wide.csv", sep = ",")

# people who didn't return for year 2

lost <- dat2[is.na(dat2$log_vload.2) & is.na(dat2$agg_ment.2),]
  #213 people lost to follow-up by year 2

dat2$hard_drugs.2 <- as.factor(dat2$hard_drugs.2)
dat2$hard_drugs.0 <- as.factor(dat2$hard_drugs.0)

Table1(lost, c("race.0", "smoke.0", "age.0", "bmi.0", "educbas.0", "hard_drugs.0"), colvar = NULL, verbose = T, incl_missing = F)
  # no startling differences from group as a whole
  # lost 27 of 66 people reporting hard drugs

############# Partial F tests ##############
# complete cases for all relevant columns
dat3 <- dat2[complete.cases(dat2[,c(7,19,26:27,29,32,57)]),]
dat3$hard_drugs.2 <- as.factor(dat3$hard_drugs.2)
dat3$hard_drugs.0 <- as.factor(dat3$hard_drugs.0)


# viral load
test1 <- lm(delta_vload ~ hard_drugs.0 + age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + log_vload.0, data = dat3)
test2 <- lm(delta_vload ~ age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + log_vload.0, data = dat3)
anova(test1, test2) # non-significant p value, hard drugs does not contribute to the outcome

# tcell
test1 <- lm(delta_tcell ~ hard_drugs.0 + age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + log_tcell.0, data = dat3)
test2 <- lm(delta_tcell ~ age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + log_tcell.0, data = dat3)
anova(test1, test2) # significant p value, hard drugs contributes to the outcome

# mental
test1 <- lm(delta_mental ~ hard_drugs.0 + age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + agg_ment.0, data = dat3)
test2 <- lm(delta_mental ~ age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + agg_ment.0, data = dat3)
anova(test1, test2) # non-significant p value, hard drugs does not contribute to the outcome

# physical
test1 <- lm(delta_phys ~ hard_drugs.0 + age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + agg_phys.0, data = dat3)
test2 <- lm(delta_phys ~ age.0 + bmi.0 + adh.2 +
              race.0 + educbas.0 + smoke.0 + agg_phys.0, data = dat3)
anova(test1, test2) # significant p value, hard drugs does contribute to the outcome


################# plotting distributions for priors #################
hist(dat2$delta_vload)
mean(dat2$delta_vload, na.rm = T) # 2.71
var(dat2$delta_vload, na.rm = T) # 1.52
  # variance is not concerning. Could keep prior variance at 1000 without much risk

hist(dat2$delta_tcell)
mean(dat2$delta_tcell, na.rm = T) # .182
  # so close to 0, it should be a sufficient prior mean
var(dat2$delta_tcell, na.rm = T) # .06
  # so small a variance of 1000 should be fine

hist(dat2$delta_mental)
mean(dat2$delta_mental, na.rm = T) # 2.21
var(dat2$delta_mental, na.rm = T) # 143.82
  # close to zero with wide variance: variance of 1000 in prior maybe not sufficient
  # close enough to 0 though that I'll leave the prior mean at 0

hist(dat2$delta_phys)
mean(dat2$delta_phys, na.rm = T) # -1.63
var(dat2$delta_phys, na.rm = T) # 71.35
  # close to 0 with narrow(er) variance. will still increase prior variance
  # will leave prior mean at 0



