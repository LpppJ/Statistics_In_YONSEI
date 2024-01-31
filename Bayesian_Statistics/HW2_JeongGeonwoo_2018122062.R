# STA3105 Bayesian Statistics
# HW2
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
dev.off()
set.seed(2018122062)

library(coda)
library(truncnorm)

##### Q1 #####
#### (a) ####
mu = 8
sigma = 1
a = 9
b = 10

truncated_normal_pdf <- function(x) {
  dnorm_part <- dnorm((x-mu)/sigma, 0, 1)
  pnorm_part <- (pnorm((b-mu)/sigma, 0, 1) - pnorm((a-mu)/sigma, 0, 1))
  indicator_part <- 1*(a <= x & x <= b)
  z_value <- (1/sigma) * dnorm_part / pnorm_part * indicator_part
  return (z_value)
}

x <- seq(9, 10, 0.001)
plot(x, truncated_normal_pdf(x), type='l', cex=1.5)

x_current = rnorm(1, mean=mu, sd=sigma)
n_samples = 10000
x_samples = rep(NA, n_samples)
acceptance = 0

while (x_current < 9 | 10 < x_current){
  x_current = rnorm(1, mean=mu, sd=sigma)
}
for (i in 1:n_samples){
  x_proposal = rnorm(1, mean=x_current, sd=1)
  U <- runif(1, 0, 1)
  alpha <- truncated_normal_pdf(x_proposal) / truncated_normal_pdf(x_current)
  if (U < alpha){
    x_new <- x_proposal
    acceptance <- acceptance + 1
  } else {
    x_new <- x_current
  }
  x_samples[i] <- x_new
  x_current <- x_new
}

x_samples
acceptance / n_samples

ts.plot(x_samples, main="TS plot of x_samples")
hist(x_samples, freq = FALSE)
lines(seq(9, 10, 0.001), truncated_normal_pdf(seq(9, 10, 0.001)), type='l')

#### (b) ####

qtruncnorm(0.25, a = 9, b = 10, mean = 8, sd = 1)
sort(rtruncnorm(10000, a = 9, b = 10, mean = 8, sd = 1))[25/100 * 10000]
sort(x_samples)[25/100 * length(x_samples)]

qtruncnorm(0.50, a = 9, b = 10, mean = 8, sd = 1)
sort(rtruncnorm(10000, a = 9, b = 10, mean = 8, sd = 1))[50/100 * 10000]
sort(x_samples)[50/100 * length(x_samples)]

qtruncnorm(0.75, a = 9, b = 10, mean = 8, sd = 1)
sort(rtruncnorm(10000, a = 9, b = 10, mean = 8, sd = 1))[75/100 * 10000]
sort(x_samples)[75/100 * length(x_samples)]


##### Q2 #####
#### (a) ####

rm(list=ls())
dev.off()
set.seed(2018122062)
library(coda)
library(truncnorm)

load(file="/Users/gunwoojung/Desktop/R_wd/Bayesian_Stat/hw02.RData")
n <- length(x)
hist(x, breaks=10)

a = 4
b = 10


log_posterior <- function(mu, sigma){
  pnorm_part <- pnorm((b-mu)/sigma) - pnorm((a-mu)/sigma)
  dnorm_part <- - (1/sigma^2) * (sum(x^2)/2 -sum(x)*mu + (n/2)*mu^2)
  prior_part <- - (mu^2)/20
  z_value <- -n*log(sigma) -n*log(pnorm_part) + dnorm_part + prior_part
  return (z_value)
}

n_samples = 10000
mu_samples = rep(NA, n_samples)
sigma_samples = rep(NA, n_samples)

mu_proposal_sd <- 5
sigma_proposal_sd <- 2

acceptance_mu <- 0
acceptance_sigma <- 0

mu_current <- rnorm(1, 0, 10)
sigma_current <- runif(1, 0, 30)

for (i in 1:n_samples){
  mu_proposal <- rnorm(1, mean=mu_current, sd=mu_proposal_sd)
  log_alpha_mu <- log_posterior(mu=mu_proposal, sigma=sigma_current) - log_posterior(mu=mu_current, sigma=sigma_current)
  if (is.nan(log_alpha_mu) | abs(log_alpha_mu)==Inf){
    mu_new <- mu_current
  } else if (log(runif(1)) < log_alpha_mu){
    mu_new <- mu_proposal
    acceptance_mu <- acceptance_mu + 1
  } else {
    mu_new <- mu_current
  }
  mu_samples[i] <- mu_new
  mu_current <- mu_new
  
  sigma_proposal <- rnorm(1, mean=sigma_current, sd=sigma_proposal_sd)
  while ((sigma_proposal < 0 | 30 < sigma_proposal)){
    sigma_proposal <- rnorm(1, mean=sigma_current, sd=sigma_proposal_sd)
  }
  log_alpha_sigma <- log_posterior(mu=mu_current, sigma=sigma_proposal) - log_posterior(mu=mu_current, sigma=sigma_current)
  if (is.nan(log_alpha_sigma) | abs(log_alpha_sigma)==Inf){
    sigma_new <- sigma_current
  } else if (log(runif(1)) < log_alpha_sigma){
    sigma_new <- sigma_proposal
    acceptance_sigma <- acceptance_sigma + 1
  } else {
    sigma_new <- sigma_current
  }
  sigma_samples[i] <- sigma_new
  sigma_current <- sigma_new
}

acceptance_mu / n_samples
acceptance_sigma / n_samples

ts.plot(mu_samples, main="TS plot of mu_samples")
ts.plot(sigma_samples, main="TS plot of sigma_samples")

hist(mu_samples, freq = FALSE, breaks = 30)
hist(sigma_samples, freq = FALSE, breaks =30)

##### Q2 #####
#### (b) ####

mean(mu_samples)
mean(sigma_samples)

hpd_interval_mu <- HPDinterval(as.mcmc(mu_samples), prob = 0.95)
hpd_interval_mu

hpd_interval_sigma <- HPDinterval(as.mcmc(sigma_samples), prob = 0.95)
hpd_interval_sigma

