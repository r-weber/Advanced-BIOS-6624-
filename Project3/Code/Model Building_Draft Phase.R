#################################
# Model Building
# Rachel Weber
# Project 3
# November 12 2018
#################################

# first get estimate of risk from all previous study periods using logistic regression
# using the predict function, predict risk for period 39
# calculate confidence intervals using bootstrapping

# library(boot)
library(tidyverse)
library(ggplot2)
library(beepr)

dat <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataRaw/vadata2.csv")
factor_vars <- c(2:4)
dat[,factor_vars] <- lapply(dat[,factor_vars] , factor)

dat <- dat[dat$proced != "2",]

table(dat$asa, dat$death30)
  # there were very few people with asa of 1 and all survived beyond 30 days. Let's combine asa 1 and 2
dat[dat$asa == "1" & !is.na(dat$asa),]$asa <- "2"

# some BMIs are unrealistic. Let's remove those above 54
# there is a bmi of 2.5, which is crazy low, but realistically not impossible. We'll keep them for now
dat[dat$bmi > 54 & !is.na(dat$bmi),]$bmi <- NA


# logistic regression
training <- dat[dat$sixmonth != "39",]
m1 <- glm(death30 ~ albumin + bmi + asa + proced, data = training)

# parametric method: repeatedly drawing from the same model, requiring more assumptions on your model
# NOT USING
#
# results1 <- matrix(nrow=44, ncol = 4)
# recent <- dat[dat$sixmonth == "39" & !is.na(dat$hospcode),]
# means1 <- matrix(nrow = 100, ncol = 3)
# 
# for(i in 1:44){
#   hospital <- recent[recent$hospcode == i,]
#   predicted <- predict(m1, newdata = hospital, type = "response")
#   for(j in 1:100){
#     s <- sample(predicted, length(predicted), replace = TRUE)
#     means1[j,1] <- round(mean(s, na.rm = T), 4)
#     means1[j,2] <- round(quantile(s, .025, na.rm =T), 4)
#     means1[j,3] <- round(quantile(s, .975, na.rm =T), 4)
#     s <- NULL
#   }
#   results1[i,2] <- round(mean(means1[,1]), 4)
#   results1[i,3] <- round(mean(means1[,2], .025, na.rm = T), 4)
#   results1[i,4] <- round(mean(means1[,3], .975, na.rm = T), 4)
#   results1[i,1] <- i
#   predicted <- NULL
# }
# colnames(results1) <- c("hospital", "mean_pred", "lower.ci", "upper.ci")
# results1 <- as.data.frame(results1)
# write.csv(results1, "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_predicted.csv")
# results1 <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_predicted.csv")



# non-parametric because your resampling data to make model every time; predictions not reliant on single sample
# runs REAL slow but runs correctly
means <- matrix(nrow = 1000, ncol = 3)
results <- matrix(nrow=44, ncol = 4)
recent <- dat[dat$sixmonth == "39" & !is.na(dat$hospcode),]

for(i in 1:44){
  hospital <- recent[recent$hospcode == i,]
  for(j in 1:1000){
    for(k in 1:44){
      d <- training[training$hospcode == k,]
      f <- sample_n(d, nrow(d), replace = T)
      if(k == 1){all <- f} else(all <- rbind(all, f))
    }
    m_loop <- glm(death30 ~ albumin + bmi + asa + proced, data = all)
    pred <- predict(m_loop, hospital, type = "response")
    means[j,1] <- round(mean(pred, na.rm = T), 4)
    s <- NULL
    all <- NULL
  }
  results[i,2] <- round(mean(means[,1]), 4)
  results[i,3] <- round(quantile(means[,1], .025, na.rm =T), 4)
  results[i,4] <- round(quantile(means[,1], .975, na.rm =T), 4)
  results[i,1] <- i
}  
colnames(results) <- c("hospital", "mean_pred", "lower.ci", "upper.ci")
results <- as.data.frame(results)
results$mean_pred <- results$mean_pred*100
results$lower.ci <- results$lower.ci*100
results$upper.ci <- results$upper.ci*100
beep(2)

write.csv(results, "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_model.csv")
results <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_model.csv")

# let's plot these CIs to see how much they truly overlap each other

ggplot(results1, aes(x = hospital, y = lower.ci)) + geom_point(color = "lightblue3") +
  geom_point(data = results1, aes(x = hospital, y = upper.ci), color = "lightblue3") + 
  geom_point(data = results, aes(x = hospital, y = lower.ci), color = "coral3") + 
  geom_point(data = results, aes(x = hospital, y = upper.ci), color = "coral3") +
  ggtitle("Lower and Upper CIs by method of Boostrap")

ggplot(results1, aes(hospital, y = mean_pred)) + geom_point(color = "lightblue3") + 
  geom_point(data = results, aes(hospital, mean_pred), color = "coral3") +
  ggtitle("Mean Predicted Rates by Method \nof Bootstrap") + ylab("Mean Predicted")
  

##################### lets recreate this without albumin in the model ###############################
means <- matrix(nrow = 150, ncol = 1)
results <- matrix(nrow = 44, ncol = 4)
recent <- dat[dat$sixmonth == "39" & !is.na(dat$hospcode),]

for(i in 1:44){
  hospital <- recent[recent$hospcode == i,]
  for(j in 1:150){
    for(k in 1:44){
      d <- training[training$hospcode == k,]
      f <- sample_n(d, nrow(d), replace = T)
      if(k == 1){all <- f} else(all <- rbind(all, f))
    }
    m_loop <- glm(death30 ~ bmi + asa + proced, data = all)
    pred <- predict(m_loop, hospital, type = "response")
    means[j,1] <- round(mean(pred, na.rm = T), 4)
    s <- NULL
    all <- NULL
  }
  results[i,2] <- round(mean(means[,1]), 4)
  results[i,3] <- round(quantile(means[,1], .025, na.rm =T), 4)
  results[i,4] <- round(quantile(means[,1], .975, na.rm =T), 4)
  results[i,1] <- i
}  
colnames(results) <- c("hospital", "mean_pred", "lower.ci", "upper.ci")
results <- as.data.frame(results)
results$mean_pred <- results$mean_pred*100
results$lower.ci <- results$lower.ci*100
results$upper.ci <- results$upper.ci*100
beep(2)

write.csv(results, "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_model_no_alb.csv")
results <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/boot_model_no_alb.csv")

########## Binomial Theorem to get variance of observed rates p. 39 ############

recent$death30 <- factor(recent$death30)
ci <- matrix(nrow = 44, ncol = 4)

for(i in 1:44){
  hospital <- recent[recent$hospcode == i,]
  n <- nrow(hospital)
  p <- nrow(hospital[hospital$death30 == 1,])/n
  q <- 1-p
  se <- sqrt(p*q/n)
  ci[i,1] <- round(p, 4)*100
  ci[i,2] <- paste(ifelse(p-1.96*se > 0, round(p-1.96*se, 4)*100, 0), "-",round(p+1.96*se, 4)*100, sep = "")
  ci[i,3] <- ifelse(p-1.96*se > 0, round(p-1.96*se, 4)*100, 0)
  ci[i,4] <- round(p+1.96*se, 4)*100
}
ci <- as.data.frame(ci)
colnames(ci) <- c("Observed", "95 CI", "lower", "upper")

write.csv(ci, "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project3/DataProcessed/binom_ci.csv")
