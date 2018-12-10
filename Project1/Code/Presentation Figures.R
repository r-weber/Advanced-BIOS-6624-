###########################
# Figures for Presentation
# Rachel Weber
# Oct 6 2018
##########################

library(ggplot2)
library(gridExtra)
library(AICcmodavg)
library(rjags)
library(mcmcse)

# read in data
dat2 <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataProcessed/hiv_wide.csv", sep = ",")
bin_col <- c(4:6, 8:13, 17, 19:22, 25:28, 30:31, 58)
dat2[, bin_col] <- lapply(dat2[, bin_col], as.factor)

# complete cases for all relevant columns
dat3 <- dat2[complete.cases(dat2[,c(7,19,26:27,29,32,57,58)]),]
dat3$hard_drugs.2 <- as.factor(dat3$hard_drugs.2)
dat3$hard_drugs.0 <- as.factor(dat3$hard_drugs.0)


# variable distributions
g1 <- ggplot(dat2, aes(x = age.0)) + geom_histogram(fill = "darkolivegreen4") +
  ggtitle("Distribution of Ages") +
  xlab("Age") + ylab("Count") + guides(fill = F) +
  theme_minimal()

g2 <- ggplot(dat2, aes(x = smoke.0, fill = smoke.0)) +
  geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = -1) +
  ggtitle("Smoking Status") +
  scale_fill_manual(values=c("orangered3", "#E69F00")) +
  xlab("Smoking Status") + ylab("Count") +
  guides(fill = F) + theme_minimal() +
  ylim(0, 600) +
  scale_x_discrete(labels = c("Non-Smoker", "Smoker"))

g3 <- ggplot(dat2, aes(x = educbas.0, fill = educbas.0)) +
  geom_bar() + geom_text(stat = 'count', aes(label = ..count..), vjust = -1) +
  ggtitle("Education Attained") +
  scale_fill_manual(values=c("orangered3", "#E69F00", "#56B4E9")) +
  xlab("Education Attained") + ylab("Count") +
  guides(fill = F) + theme_minimal() +
  ylim(0, 600) +
  scale_x_discrete(labels = c("< College", "College Degree", "> College"))

grid.arrange(g1, g2, g3, nrow = 3)


# distributions of outcomes

g4 <- ggplot(dat3, aes(x = delta_mental)) + geom_histogram(color = "#56B4E9", fill = "#56B4E9") +
  ggtitle("Change in Mental Health Score") + xlab("Score Change") +
  theme_minimal()
g5 <- ggplot(dat3, aes(x = delta_phys)) + geom_histogram(color = "#E69F00", fill = "#E69F00") +
  ggtitle("Change in Physical Health Score") + xlab("Score Change") +
  theme_minimal()
g6 <- ggplot(dat3, aes(x = delta_vload)) + geom_histogram(color = "orangered3", fill = "orangered3") +
  ggtitle("Change in log(Viral Load)") + xlab("Change in Viral Load") +
  theme_minimal()
g7 <- ggplot(dat3, aes(x = delta_tcell)) + geom_histogram(color = "darkolivegreen4", fill = "darkolivegreen4") +
  ggtitle("Change in log(T Cell Count)") + xlab("Change in T Cell Count") +
  theme_minimal()

grid.arrange(g4, g5, g6, g7, nrow = 4)


dat3 <- read.csv(file = "C:/Users/weberra/Documents/Classwork/Advanced/bios6624-r-weber/Project1/DataProcessed/hivfinal.csv", sep = ",")

g4 <- ggplot(dat3, aes(x = delta.agg.ment)) + geom_histogram(color = "#56B4E9", fill = "#56B4E9") +
  ggtitle("Change in Mental Health Score") + xlab("Score Change") +
  theme_minimal()
g5 <- ggplot(dat3, aes(x = delta.agg.phys)) + geom_histogram(color = "#E69F00", fill = "#E69F00") +
  ggtitle("Change in Physical Health Score") + xlab("Score Change") +
  theme_minimal()
g6 <- ggplot(dat3, aes(x = delta.log.viralload)) + geom_histogram(color = "orangered3", fill = "orangered3") +
  ggtitle("Change in log(Viral Load)") + xlab("Change in Viral Load") +
  theme_minimal()
g7 <- ggplot(dat3, aes(x = delta.log.leu3n)) + geom_histogram(color = "darkolivegreen4", fill = "darkolivegreen4") +
  ggtitle("Change in log(T Cell Count)") + xlab("Change in T Cell Count") +
  theme_minimal()

grid.arrange(g4, g5, g6, g7, nrow = 4)


###################### trace plots and autocorrelation #####################
# model of T-Cell with hard drugs
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
iter <- 25000
samples <- coda.samples(mod2, variable.names = c("beta", "sigma2"), n.iter = iter)

# Explore support of the parameter space
samples_dic <- dic.samples(mod2, n.iter = iter, type = "pD")
plot(samples)
# looks good. Distribution seems well explored

# model of Physical QOL with hard drugs
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

plot(samples)