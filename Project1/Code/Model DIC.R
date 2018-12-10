##########################
# DIC for Bayes Models
# Rachel Weber
# Oct 6 18
#########################

# results summary:
# viral load: model without is better...(1433 vs 1431)
# T-Cell: model with is slightly better...(-201 vs -165)
# mental health: model without hard drugs is better...(3442 vs 3441)
# physical health: model with is better (3223 vs 3227)



library(AICcmodavg)
library(rjags)
library(mcmcse)

dat2 <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataProcessed/hiv_wide.csv", sep = ",")
bin_col <- c(4:6, 8:13, 17, 19:22, 25:28, 30:31, 58)
dat2[, bin_col] <- lapply(dat2[, bin_col], as.factor)

# complete cases for all relevant columns
dat3 <- dat2[complete.cases(dat2[,c(7,19,26:27,29,32,57,58)]),]
dat3$hard_drugs.2 <- as.factor(dat3$hard_drugs.2)
dat3$hard.drogs0 <- as.factor(dat3$hard.drugs0)

# alright, run all of the models and then get DIC for everybody

dat3 <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataProcessed/hivfinal.csv", sep = ",")


################################ Viral Load ##################################
# model with hard drugs
y <- c(dat3$delta.log.viralload)
X <- model.matrix(~ hard.drugs0 + age + bmi + adh +
                    race + educ + smoke + log.viralload0, data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.000001
b <- 0.000001
m <- rep(2.71, p) # this is where the mean comes in... changing from 0 to 2.71
R <- matrix(0, p, p)
diag(R) <- 0.000001 # widening variance to account for big outliers

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)
mod1 <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

dic.samples(mod1, n.iter = 1000, thin = 1)

iter <- 25000
samples <- coda.samples(mod1, variable.names = c("beta", "sigma2"), n.iter = iter)
hpd <- function(x, alpha = 0.05){
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
}

draws <- as.matrix(samples)
out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

pander(out_mat)


# model without hard drugs
y <- c(dat3$delta.log.viralload)
X <- model.matrix(~ age + bmi + adh +
                    race + educ + smoke + log.viralload0, data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.000001
b <- 0.000001
m <- rep(2.71, p) # this is where the mean comes in... changing from 0 to 2.71
R <- matrix(0, p, p)
diag(R) <- 0.000001 # widening variance to account for big outliers

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)
mod1_nohd <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

dic.samples(mod1_nohd, n.iter = 1000, thin = 1)
# model without hard drugs is *slightly* worse.

################################## TCell count ##################################
# model with hard drugs
y <- c(dat3$delta.log.leu3n)

X <- model.matrix(~ hard.drugs0 + age + bmi + adh +
                    race + educ + smoke + log.leu3n0, data = dat3)
N <- nrow(X)
p <- ncol(X)

a <- 0.0001
b <- 0.0001
m <- rep(0, p) # this is where the mean is specified, changed from 0 to .182
R <- matrix(0, p, p)
diag(R) <- 0.0001

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)

mod2 <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

dic.samples(mod2, n.iter = 1000, thin = 1)

# Explore support of the parameter space
samples_dic <- dic.samples(mod2, n.iter = iter, type = "pD")
plot(samples)
# looks good. Distribution seems well explored

# model without hard drugs
y <- c(dat3$delta.log.leu3n)
X <- model.matrix(~ age + bmi + adh +
                    race + educ + smoke + log.leu3n0, data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.0001
b <- 0.0001
m <- rep(0, p) # this is where the mean is specified, changed from 0 to .182
R <- matrix(0, p, p)
diag(R) <- 0.0001

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)
mod2_nhd <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

dic.samples(mod2_nhd, n.iter = 1000, thin = 1)

######################################## mental health ###########################################


# model with hard drugs
y <- c(dat3$delta.agg.ment)
X <- model.matrix(~ hard.drugs0 + age + bmi + adh + race + 
                        educ + smoke + agg.ment0 , data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.000001
b <- 0.000001
m <- rep(0, p)
R <- matrix(0, p, p)
diag(R) <- 0.000001 # widening variance for good measure

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)
mod3 <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

iter <- 25000
samples <- coda.samples(mod3, variable.names = c("beta", "sigma2"), n.iter = iter)
hpd <- function(x, alpha = 0.05){
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
}

draws <- as.matrix(samples)
out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

pander(out_mat)

dic.samples(mod3, n.iter = 1000, thin = 1)

# Explore support of the parameter space
samples_dic <- dic.samples(mod3, n.iter = iter, type = "pD")
plot(samples)
# looks good. Distribution seems well explored

# model without hard drugs
y <- c(dat3$delta.agg.ment)
X <- model.matrix(~ age + bmi + adh +
                    race + educ + smoke + agg.ment0, data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.000001
b <- 0.000001
m <- rep(0, p)
R <- matrix(0, p, p)
diag(R) <- 0.000001 # widening variance for good measure

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)
mod3_nohd <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

dic.samples(mod3_nohd, n.iter = 1000, thin = 1)

################################# physical health ################################
# model with hard drugs
y <- c(dat3$delta.agg.phys)
X <- model.matrix(~ hard.drugs0 + age + bmi + adh + race + 
                    educ + smoke + agg.phys0 , data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.000001
b <- 0.000001
m <- rep(0, p)
R <- matrix(0, p, p)
diag(R) <- 0.000001 # widening variance for good measure

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)

mod4 <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

samples <- coda.samples(mod4, variable.names = c("beta", "sigma2"), n.iter = iter)

draws <- as.matrix(samples)
out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

pander(out_mat)

dic.samples(mod4, n.iter = 1000, thin = 1)

# explore support of the parameter space
samples_dic <- dic.samples(mod, n.iter = iter, type = "pD")
plot(samples)
# looks good. Distribution seems well explored

# model without hard drugs
y <- c(dat3$delta.agg.phys)
X <- model.matrix(~ age + bmi + adh +
                    race + educ + smoke + agg_phys.0, data = dat3)
N <- nrow(X)
p <- ncol(X)
a <- 0.000001
b <- 0.000001
m <- rep(0, p)
R <- matrix(0, p, p)
diag(R) <- 0.000001 # widening variance for good measure

jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)

mod4_nohd <- jags.model("C:/Users/weberra/Documents/Classwork/Advanced/InClass/Bayesian/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

dic.samples(mod4_nohd, n.iter = 1000, thin = 1)