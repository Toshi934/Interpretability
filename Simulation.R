#######Simulation Code#######
###Last Update: 2024/08/21###
rm(list = ls(all.names = TRUE))

#########################################################################
###Install packages
library(tidyverse)
library(bcf)
library(lmtest)
library(sandwich)
library(rpart)


#########################################################################
###Test data generation

#Sample size
N <- 10000

#Treatment
set.seed(1)
treatment <- rbinom(N, 1, 0.5)

#Covariates
age <- rnorm(N, 40, 10) 
systolic.blood.pressure <- rnorm(N, 120, 10) 
hba1c <- rnorm(N, 5, 1) 
eGFR <- rnorm(N, 60, 20) 
eGFR <- pmin(pmax(eGFR, 10), 120)
medication <- rbinom(N, 1, 0.5)

#Parameters for HTE
b <- (0.1 + 0.05*hba1c + 0.05*eGFR + 0.3*medication) / 20

#Outcome
regY1 <- -0.1 + b*treatment + 0.0005*(systolic.blood.pressure - mean(systolic.blood.pressure)) + 0.01*hba1c + 0.001*(eGFR- mean(eGFR)) + 0.01*medication
regY1 <- pmin(pmax(regY1, 0), 1)
outcome <- rbinom(N, 1, regY1)

#Propensity score
pi <- rep(0.5, N)

#Create a data frame
data <- data.frame(treatment, outcome, age, systolic.blood.pressure, hba1c,  eGFR, medication, pi)

#########################################################################
###Step1: CATE estimation

###Run Bayesian causal forest###
Y <-  data$outcome
W <-  data$treatment
PS <- data$pi
X <- data[, c("age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")]
X <- as.matrix(X)

bcf.model <- bcf(y = Y, #Outcome
                 z = W, #Treatment
                 x_control = X, #Covariate for mu (prognostic score)
                 x_moderate = X, #Covariate for tau (treatment effect)
                 pihat = PS, #Propensity score
                 nburn = 300, #The number of iterations for burn-in
                 nsim = 300, #The number of iterations used after burn-in 
                 ntree_control = 200, #The number of trees used for mu estimation 
                 ntree_moderate = 200, #The number of trees used for tau estimation s
                 random_seed = 1) #Random seed 

#Obtain estimated CATE
tau <- bcf.model$tau
data$tau.hat <- colMeans(tau)


###Test calibration performance of Bayesian causal forest###
#Convert tau from binary to a percentage scale
data$tau.hat <- data$tau.hat * 100

#Classify sample by CATE quintile
data$ranking <- ntile(data$tau.hat, 5) %>% as.factor

#Estimate ATE for each CATE subgroup (linear regression)
ols.ate <- lm(outcome ~ 0 + ranking + ranking:treatment, data = data)
ols.ate <- coeftest(ols.ate, vcov=vcovHC(ols.ate, type = "HC2"))

#Save results of linear regression
interact <- which(grepl (":", rownames(ols.ate)))
ols.ate <- data.frame("OLS", paste0("Q", seq(5)), ols.ate[interact, 1:2])
rownames (ols.ate) <- NULL
colnames (ols.ate) <- c("method", "ranking", "estimate", "std.err")
ols.ate$ranking <- factor(ols.ate$ranking, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

#Plot estimated ATEs by CATE quintile
calibration.plot <- ggplot(ols.ate) +
  aes (x = ranking, y = estimate, group = method, color = method) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = estimate -2 * std.err, ymax = estimate +2 * std.err), width = .2, position = position_dodge (0.2)) +
  ylab("ATE (percentage point)") + xlab("CATE Ranking") + 
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) 

print(calibration.plot)


#########################################################################
###Step2: subgroup identification by CARR

###Create interpretable covariates###
#We convert all continuous variables into categories 
data.cart <- data
data.cart <- data.cart %>% mutate(age = case_when(
  age < 40 ~ 1, #Low
  age >= 40 & age < 60 ~ 2, #Middle
  age >= 60 ~ 3 #High
  ),
  systolic.blood.pressure = case_when(
    systolic.blood.pressure < 120 ~ 1, 
    systolic.blood.pressure >= 120 & systolic.blood.pressure < 130 ~ 2, 
    systolic.blood.pressure >= 130 & systolic.blood.pressure < 140 ~ 3, 
    systolic.blood.pressure >= 140 ~ 4
  ),
  hba1c = case_when(
    hba1c < 5.7 ~ 1, 
    hba1c >= 5.7 & hba1c < 6.5 ~ 2, 
    hba1c >= 6.5 ~ 3 
  ),
  eGFR = case_when(
    eGFR >= 90 ~ 1, 
    eGFR < 90 & eGFR >= 60 ~ 2, 
    eGFR < 60 ~ 3
  ))


###Run CART to identify subgroups###
cart.model <- rpart(tau.hat ~ age + systolic.blood.pressure + hba1c + eGFR + medication, 
                    data = data.cart, 
                    maxdepth = 2)

#Plot the identified subgroups
par(xpd = NA)
plot(cart.model, uniform = TRUE)
text(cart.model, digits = 3, use.n = TRUE)

#Result indicates samples with 
#1) eGFR <60 & not taking medication have the lowest level of CATE
#2) eGFR ≥90  have the highest level of CATE


