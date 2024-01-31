# STA3105 Bayesian Statistics
# HW2
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)
library(MASS)
library(ggplot2)
library(coda)
##### Q1 #####
X <- mvrnorm(n = 1000, mu = rep(0, 4), Sigma = diag(4))
beta_true <- c(0.5, -0.5, 0, 1)

# Logistic function
mu <- function(X, beta) {
  exp_Xbeta <- exp(X %*% beta)
  return(exp_Xbeta / (1 + exp_Xbeta))
}

mu_X <- mu(X, beta_true)
Y <- rbinom(1000, 10, mu_X)

# hist(mu_X, breaks = 15)
# hist(Y, breaks = -1:10)

##### Q2 #####

log_likelihood <- function(Y, X, beta) {
  p <- mu(X, beta)
  return(sum(Y * log(p) + (10 - Y) * log(1 - p)))
}

log_prior <- function(beta) {
  return(-beta^2 / 20)
}

# Run the Metropolis-Hastings Algorithm
iterations <- 10000
proposal_sd <- rep(0.1,4)
beta_samples <- matrix(NA, nrow = iterations, ncol = ncol(X))
# initialization
current_beta <- rep(0,4)
beta_samples[1, ] <- current_beta
# the number of acceptance
acceptance = rep(0,4)

for(i in 2:iterations) {
  for (j in 1:ncol(X)){
    # Propose a beta_j
    proposed_beta <- current_beta
    proposed_beta[j] <- rnorm(1, current_beta[j], proposal_sd[j])
    # Calculate acceptance probability
    proposed_log_posterior <- log_likelihood(Y, X, proposed_beta) + log_prior(proposed_beta[j])
    current_log_posterior <- log_likelihood(Y, X, current_beta) + log_prior(current_beta[j])
    log_alpha <- proposed_log_posterior - current_log_posterior
    
    # Accept or reject the new beta
    if (is.nan(log_alpha) | abs(log_alpha)==Inf){
    # if at least one posterior is nan or inf, -inf
    # reject
      
    } else if(log(runif(1)) < log_alpha) {
      # the condition of accept
      current_beta[j] <- proposed_beta[j]
      acceptance[j] = acceptance[j] + 1
    }
    beta_samples[i, ] <- current_beta
  }
}

# Burn-in (After Checking TS plot, then I modified!)
burn_in <- 100
beta_samples <- beta_samples[seq(burn_in, nrow(beta_samples)), ]

# ts plots
# par(mfrow=c(2,2))
ts.plot(beta_samples[,1], main="TS plot of beta1", xlab="iteration", ylab="beta1")
ts.plot(beta_samples[,2], main="TS plot of beta2", xlab="iteration", ylab="beta2")
ts.plot(beta_samples[,3], main="TS plot of beta3", xlab="iteration", ylab="beta3")
ts.plot(beta_samples[,4], main="TS plot of beta4", xlab="iteration", ylab="beta4")

# density plots
# par(mfrow=c(2,2))
plot(density(beta_samples[,1]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples[,2]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples[,3]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples[,4]), main="Density Plot of beta4", xlab="beta", ylab="Density", lwd=2)

# HPD intervals
cat("95% HPD intervals of beta1:", HPDinterval(as.mcmc(beta_samples[,1]), prob = 0.95), "\n",
    "95% HPD intervals of beta2:", HPDinterval(as.mcmc(beta_samples[,2]), prob = 0.95), "\n",
    "95% HPD intervals of beta3:", HPDinterval(as.mcmc(beta_samples[,3]), prob = 0.95), "\n",
    "95% HPD intervals of beta4:", HPDinterval(as.mcmc(beta_samples[,4]), prob = 0.95))

# Posterior mean
cat("posterior mean of beta1:", mean(beta_samples[,1]))
cat("posterior mean of beta2:", mean(beta_samples[,2]))
cat("posterior mean of beta3:", mean(beta_samples[,3]))
cat("posterior mean of beta4:", mean(beta_samples[,4]))

# Acceptance rate
cat("acceptance probability of beta1:", length(unique(beta_samples[,1]))/iterations)
cat("acceptance probability of beta2:", length(unique(beta_samples[,2]))/iterations)
cat("acceptance probability of beta3:", length(unique(beta_samples[,3]))/iterations)
cat("acceptance probability of beta4:", length(unique(beta_samples[,4]))/iterations)

# Thinning (if we need)
# thinning = 2
# beta_samples <- beta_samples[seq(1, nrow(beta_samples), by=thinning), ]

# Effective sample size (after thinning)
cat("effective sample size of beta1:", effectiveSize(beta_samples[,1]))
cat("effective sample size of beta2:", effectiveSize(beta_samples[,2]))
cat("effective sample size of beta3:", effectiveSize(beta_samples[,3]))
cat("effective sample size of beta4:", effectiveSize(beta_samples[,4]))

##### Q3 #####
graphics.off()
library(nimble)
model_string <- nimbleCode({
  for(i in 1:length(Y)) {
    logit(p[i]) <- inprod(X[i, 1:4], beta[1:4])
    Y[i] ~ dbin(p[i], 10)
  }
  for(j in 1:4) {
    beta[j] ~ dnorm(0, sd = sqrt(10))
  }
})

model <- nimbleModel(
  code = model_string,
  data = list(Y = Y),
  constants = list(X = X),
  inits = list(beta = rep(0,4)),
  name = "LogisticRegression"
)

compiled_model <- compileNimble(model)
mcmc <- buildMCMC(compiled_model)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
results <- runMCMC(compiled_mcmc, niter = iterations, nburnin = burn_in, thin = 1)
beta_samples_nimble <- as.matrix(results, 'beta')

# par(mfrow=c(2,2))
ts.plot(beta_samples_nimble[,1], main="TS plot of beta1", xlab="iteration", ylab="beta1")
ts.plot(beta_samples_nimble[,2], main="TS plot of beta2", xlab="iteration", ylab="beta2")
ts.plot(beta_samples_nimble[,3], main="TS plot of beta3", xlab="iteration", ylab="beta3")
ts.plot(beta_samples_nimble[,4], main="TS plot of beta4", xlab="iteration", ylab="beta4")

# par(mfrow=c(2,2))
plot(density(beta_samples_nimble[,1]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples_nimble[,2]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples_nimble[,3]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples_nimble[,4]), main="Density Plot of beta4", xlab="beta", ylab="Density", lwd=2)

hpd_interval_beta1 <- HPDinterval(as.mcmc(beta_samples_nimble[,1]), prob = 0.95)
cat("95% HPD intervals of beta1:", hpd_interval_beta1)
hpd_interval_beta2 <- HPDinterval(as.mcmc(beta_samples_nimble[,2]), prob = 0.95)
cat("95% HPD intervals of beta2:", hpd_interval_beta2)
hpd_interval_beta3 <- HPDinterval(as.mcmc(beta_samples_nimble[,3]), prob = 0.95)
cat("95% HPD intervals of beta3:", hpd_interval_beta3)
hpd_interval_beta4 <- HPDinterval(as.mcmc(beta_samples_nimble[,4]), prob = 0.95)
cat("95% HPD intervals of beta4:", hpd_interval_beta4)

cat("posterior mean of beta1:", mean(beta_samples_nimble[,1]))
cat("posterior mean of beta2:", mean(beta_samples_nimble[,2]))
cat("posterior mean of beta3:", mean(beta_samples_nimble[,3]))
cat("posterior mean of beta4:", mean(beta_samples_nimble[,4]))

cat("acceptance probability of beta1:", length(unique(beta_samples_nimble[,1]))/iterations)
cat("acceptance probability of beta2:", length(unique(beta_samples_nimble[,2]))/iterations)
cat("acceptance probability of beta3:", length(unique(beta_samples_nimble[,3]))/iterations)
cat("acceptance probability of beta4:", length(unique(beta_samples_nimble[,4]))/iterations)

cat("effective sample size of beta1:", effectiveSize(beta_samples_nimble[,1]))
cat("effective sample size of beta2:", effectiveSize(beta_samples_nimble[,2]))
cat("effective sample size of beta3:", effectiveSize(beta_samples_nimble[,3]))
cat("effective sample size of beta4:", effectiveSize(beta_samples_nimble[,4]))

##### Q4 #####
graphics.off()
library(adaptMCMC)
init.pars = c(beta1=0, beta2=0, beta3=0, beta4=0)

# Define log posterior
log_posterior <- function(pars) {
  with(as.list(pars),{
    beta = pars
    logit_p <- X %*% beta
    p <- exp(logit_p) / (1 + exp(logit_p))
    log_prior <- sum(dnorm(beta, 0, sqrt(10), log = TRUE))
    log_lik <- sum(dbinom(Y, 20, p, log = TRUE))
    return (log_prior + log_lik)
  })
}

iterations = 13000
mcmc_results <- MCMC(p=log_posterior, n=iterations, init=rep(0,4), scale=rep(10,4), adapt=TRUE, acc.rate=0.3)
beta_samples_adaptMCMC <- mcmc_results$samples

# Burn-in (After Checking TS plot, then I modified!)
burn_in <- 3000
beta_samples_adaptMCMC <- beta_samples_adaptMCMC[seq(burn_in, nrow(beta_samples_adaptMCMC)), ]

# ts plots
# par(mfrow=c(2,2))
ts.plot(beta_samples_adaptMCMC[,1], main="TS plot of beta1", xlab="iteration", ylab="beta1")
ts.plot(beta_samples_adaptMCMC[,2], main="TS plot of beta2", xlab="iteration", ylab="beta2")
ts.plot(beta_samples_adaptMCMC[,3], main="TS plot of beta3", xlab="iteration", ylab="beta3")
ts.plot(beta_samples_adaptMCMC[,4], main="TS plot of beta4", xlab="iteration", ylab="beta4")

# par(mfrow=c(2,2))
plot(density(beta_samples_adaptMCMC[,1]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples_adaptMCMC[,2]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples_adaptMCMC[,3]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samples_adaptMCMC[,4]), main="Density Plot of beta4", xlab="beta", ylab="Density", lwd=2)

hpd_interval_beta1 <- HPDinterval(as.mcmc(beta_samples_adaptMCMC[,1]), prob = 0.95)
cat("95% HPD intervals of beta1:", hpd_interval_beta1)
hpd_interval_beta2 <- HPDinterval(as.mcmc(beta_samples_adaptMCMC[,2]), prob = 0.95)
cat("95% HPD intervals of beta2:", hpd_interval_beta2)
hpd_interval_beta3 <- HPDinterval(as.mcmc(beta_samples_adaptMCMC[,3]), prob = 0.95)
cat("95% HPD intervals of beta3:", hpd_interval_beta3)
hpd_interval_beta4 <- HPDinterval(as.mcmc(beta_samples_adaptMCMC[,4]), prob = 0.95)
cat("95% HPD intervals of beta4:", hpd_interval_beta4)

cat("posterior mean of beta1:", mean(beta_samples_adaptMCMC[,1]))
cat("posterior mean of beta2:", mean(beta_samples_adaptMCMC[,2]))
cat("posterior mean of beta3:", mean(beta_samples_adaptMCMC[,3]))
cat("posterior mean of beta4:", mean(beta_samples_adaptMCMC[,4]))

cat("acceptance probability of beta1:", length(unique(beta_samples_adaptMCMC[,1]))/iterations)
cat("acceptance probability of beta2:", length(unique(beta_samples_adaptMCMC[,2]))/iterations)
cat("acceptance probability of beta3:", length(unique(beta_samples_adaptMCMC[,3]))/iterations)
cat("acceptance probability of beta4:", length(unique(beta_samples_adaptMCMC[,4]))/iterations)

# Thinning (if we need)
# thinning = 2
# beta_samples_adaptMCMC <- beta_samples_adaptMCMC[seq(1, nrow(beta_samples_adaptMCMC), by=thinning), ]

cat("effective sample size of beta1:", effectiveSize(beta_samples_adaptMCMC[,1]))
cat("effective sample size of beta2:", effectiveSize(beta_samples_adaptMCMC[,2]))
cat("effective sample size of beta3:", effectiveSize(beta_samples_adaptMCMC[,3]))
cat("effective sample size of beta4:", effectiveSize(beta_samples_adaptMCMC[,4]))

##### Q5 #####
dev.off()

for (i in 1:4){
  x_lim = range(c(
    density(beta_samples[, i])$x,
    density(beta_samples_nimble[, i])$x, 
    density(beta_samples_adaptMCMC[, i])$x
  ))
  y_lim = range(c(
    density(beta_samples[, i])$y, 
    density(beta_samples_nimble[, i])$y, 
    density(beta_samples_adaptMCMC[, i])$y
  ))
  # Initial plot settings
  plot(
    density(beta_samples[, i]), 
    main = paste("Density Plot of beta", i),
    xlab = "beta1", 
    ylab = "Density", 
    lwd = 2,
    col = "red", 
    xlim = x_lim, 
    ylim = y_lim
  )
  
  lines(density(beta_samples_nimble[, i]), lwd = 2, col = "blue")
  lines(density(beta_samples_adaptMCMC[, i]), lwd = 2, col = "green")
  abline(v = beta_true[i], lwd = 2, col = "black")
  
  legend_location_x = c(0.3, -0.45, 0.01, 0.6)
  legend(x = legend_location_x[i], y = max(y_lim)+1,  # x와 y 좌표를 조정하여 위치 지정
    legend = c("beta_samples_no_packages", "beta_samples_nimble", "beta_samples_adaptMCMC","beta_true"),
    col = c("red", "blue", "green","black"),
    lwd = 2, cex = 0.7, box.lwd = 0, bg = "transparent",
  )
}