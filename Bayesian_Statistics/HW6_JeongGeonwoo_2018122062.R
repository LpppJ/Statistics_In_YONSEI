# STA3105 Bayesian Statistics
# HW6
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)
library(adaptMCMC); library(coda); library(invgamma); library(mvtnorm)

##### Q1 #####

load("/Users/gunwoojung/Desktop/R_wd/Bayesian_Stat/hw6.RData")
X.train <- X[1:70]; Y.train <- Y[1:70]
X.test <- X[71:100]; Y.test <- Y[71:100]

model1.log.posterior <- function(pars){
  with(as.list(pars),{
    beta <- pars[1:2]; sigma2 <- pars[3]
    log.prior.beta <- sum(dnorm(beta, 0, 10, log=TRUE))
    log.prior.sigma2 <- log(dinvgamma(sigma2, 3, 3))
    Y.hat <- beta[1] + beta[2] * X.train
    log.lik <- sum(dnorm(Y.train, mean=Y.hat, sd=sqrt(sigma2), log=TRUE))
    log.posterior <- log.prior.beta + log.prior.sigma2 + log.lik
    return (log.posterior)
  })
}

model2.log.posterior <- function(pars){
  with(as.list(pars),{
    beta <- pars[1:2]; sigma2 <- pars[3]
    log.prior.beta <- sum(dnorm(beta, 0, 10, log=TRUE))
    log.prior.sigma2 <- log(dinvgamma(sigma2, 3, 3))
    Y.hat <- beta[1] + beta[2] * X.train^3
    log.lik <- sum(dnorm(Y.train, mean=Y.hat, sd=sqrt(sigma2), log=TRUE))
    log.posterior <- log.prior.beta + log.prior.sigma2 + log.lik
    return (log.posterior)
  })
}

iter <- 10000; init.pars <- c(beta0=1, beta1=1, sigma2=1)
model1.mcmc <- MCMC(p=model1.log.posterior, n=iter, init=init.pars, adapt=TRUE, acc.rate=0.4)
model2.mcmc <- MCMC(p=model2.log.posterior, n=iter, init=init.pars, adapt=TRUE, acc.rate=0.4)

# trace plot (for checking burn-in)
par(mfrow=c(1,3))
ts.plot(model1.mcmc$samples[,1], main="Trace plot of beta0 (M1)")
ts.plot(model1.mcmc$samples[,2], main="Trace plot of beta1 (M1)")
ts.plot(model1.mcmc$samples[,3], main="Trace plot of sigma2 (M1)")

par(mfrow=c(1,3))
ts.plot(model2.mcmc$samples[,1], main="Trace plot of beta0 (M2)")
ts.plot(model2.mcmc$samples[,2], main="Trace plot of beta1 (M2)")
ts.plot(model2.mcmc$samples[,3], main="Trace plot of sigma2 (M2)")

# burn in
burn.in.size = 1000
model1.mcmc$samples <- model1.mcmc$samples[(burn.in.size+1):iter,]
model2.mcmc$samples <- model2.mcmc$samples[(burn.in.size+1):iter,]

# trace plot (after burn-in)
par(mfrow=c(1,3))
ts.plot(model1.mcmc$samples[,1], main="Trace plot of beta0 (M1)")
ts.plot(model1.mcmc$samples[,2], main="Trace plot of beta1 (M1)")
ts.plot(model1.mcmc$samples[,3], main="Trace plot of sigma2 (M1)")

par(mfrow=c(1,3))
ts.plot(model2.mcmc$samples[,1], main="Trace plot of beta0 (M2)")
ts.plot(model2.mcmc$samples[,2], main="Trace plot of beta1 (M2)")
ts.plot(model2.mcmc$samples[,3], main="Trace plot of sigma2 (M2)")

# density plot
par(mfrow=c(1,3))
plot(density(model1.mcmc$samples[,1]), main="Trace plot of beta0 (M1)")
plot(density(model1.mcmc$samples[,2]), main="Trace plot of beta1 (M1)")
plot(density(model1.mcmc$samples[,3]), main="Trace plot of sigma2 (M1)")

par(mfrow=c(1,3))
plot(density(model2.mcmc$samples[,1]), main="Trace plot of beta0 (M2)")
plot(density(model2.mcmc$samples[,2]), main="Trace plot of beta1 (M2)")
plot(density(model2.mcmc$samples[,3]), main="Trace plot of sigma2 (M2)")

# 95% HPD intervals
cat("95% HPD intervals of beta0 (M1):", HPDinterval(as.mcmc(model1.mcmc$samples[,1]), prob=0.95), "\n",
    "95% HPD intervals of beta1 (M1):", HPDinterval(as.mcmc(model1.mcmc$samples[,2]), prob=0.95), "\n",
    "95% HPD intervals of sigma2 (M1):", HPDinterval(as.mcmc(model1.mcmc$samples[,3]), prob=0.95))

cat("95% HPD intervals of beta0 (M2):", HPDinterval(as.mcmc(model2.mcmc$samples[,1]), prob=0.95), "\n",
    "95% HPD intervals of beta1 (M2):", HPDinterval(as.mcmc(model2.mcmc$samples[,2]), prob=0.95), "\n",
    "95% HPD intervals of sigma2 (M2):", HPDinterval(as.mcmc(model2.mcmc$samples[,3]), prob=0.95))

# posterior mean
cat("posterior mean of beta0 (M1):", mean(model1.mcmc$samples[,1]), "\n",
    "posterior mean of beta1 (M1):", mean(model1.mcmc$samples[,2]), "\n",
    "posterior mean of sigma2 (M1):", mean(model1.mcmc$samples[,3]))

cat("posterior mean of beta0 (M2):", mean(model2.mcmc$samples[,1]), "\n",
    "posterior mean of beta1 (M2):", mean(model2.mcmc$samples[,2]), "\n",
    "posterior mean of sigma2 (M2):", mean(model2.mcmc$samples[,3]))

# acceptance probability
cat("acceptance probability of beta0 (M1):", length(unique(model1.mcmc$samples[,1]))/iter, "\n",
    "acceptance probability of beta1 (M1):", length(unique(model1.mcmc$samples[,2]))/iter, "\n",
    "acceptance probability of sigma2 (M1):", length(unique(model1.mcmc$samples[,3]))/iter)

cat("acceptance probability of beta0 (M2):", length(unique(model2.mcmc$samples[,1]))/iter, "\n",
    "acceptance probability of beta1 (M2):", length(unique(model2.mcmc$samples[,2]))/iter, "\n",
    "acceptance probability of sigma2 (M2):", length(unique(model2.mcmc$samples[,3]))/iter)

# Effective sample size
cat("effective sample size of beta0 (M1):", effectiveSize(model1.mcmc$samples[,1]), "\n",
    "effective sample size of beta1 (M1):", effectiveSize(model1.mcmc$samples[,2]), "\n",
    "effective sample size of sigma2 (M1):", effectiveSize(model1.mcmc$samples[,3]))

cat("effective sample size of beta0 (M2):", effectiveSize(model2.mcmc$samples[,1]), "\n",
    "effective sample size of beta1 (M2):", effectiveSize(model2.mcmc$samples[,2]), "\n",
    "effective sample size of sigma2 (M2):", effectiveSize(model2.mcmc$samples[,3]))

##### Q2 #####
M1.beta0.mean <- mean(model1.mcmc$samples[,1])
M1.beta1.mean <- mean(model1.mcmc$samples[,2])
M1.sigma2.mean <- mean(model1.mcmc$samples[,3])

M2.beta0.mean <- mean(model2.mcmc$samples[,1])
M2.beta1.mean <- mean(model2.mcmc$samples[,2])
M2.sigma2.mean <- mean(model2.mcmc$samples[,3])


M1.Y.hat <- M1.beta0.mean + M1.beta1.mean * X.train
M1.BIC <- -2 * dmvnorm(Y.train, M1.Y.hat, M1.sigma2.mean * diag(length(Y.train)), log = TRUE
                   ) + ncol(model1.mcmc$samples) * log(length(X.train))

M2.Y.hat <- M2.beta0.mean + M2.beta1.mean * X.train^3
M2.BIC <- -2 * dmvnorm(Y.train, M2.Y.hat, M2.sigma2.mean * diag(length(Y.train)), log = TRUE
                   ) + ncol(model2.mcmc$samples) * log(length(X.train))

cat("BIC of M1:", M1.BIC, "||", "BIC of M2:", M2.BIC)

##### Q3 #####
# model prior
M1.prior <- 0.5; M2.prior <- 0.5

# model evidence
M1.evd <- exp(-0.5 * M1.BIC); M2.evd <- exp(-0.5 * M2.BIC)

# Posterior model probability (PMP)
Sum.weight <- M1.evd * M1.prior + M2.evd * M2.prior
M1.weight <- M1.evd * M1.prior / Sum.weight
M2.weight <- M2.evd * M2.prior / Sum.weight

# Posterior predictive distribution Test[1:10]
# I defined a Empty Space "PPD.Test"
# Shape : (iter - burn.in.size) * (length(X.test)) = 9000 * 30
# Ex. Posterior predictive distribution of Test[k] (k=1,...,30)
#     will be saved k-th column of PPT.Test !
PPD.Test <- matrix( rep(NA, (iter - burn.in.size) * length(X.test)) ,
                    ncol=length(X.test))

for (k in 1:length(X.test)){
  M1.PPD.Test.k <- model1.mcmc$samples[,1] + model1.mcmc$samples[,2] * X.test[k] + rnorm(9000, 0, sd=sqrt(model1.mcmc$samples[,3][k]))
  M2.PPD.Test.k <- model2.mcmc$samples[,1] + model2.mcmc$samples[,2] * X.test[k]^3 + rnorm(9000, 0, sd=sqrt(model2.mcmc$samples[,3][k]))
  PPD.Test[,k] = M1.weight * M1.PPD.Test.k + M2.weight * M2.PPD.Test.k}

for (i in seq(1, 9, 2)) {
  par(mfrow = c(1, 2))
  plot(density(PPD.Test[, i]), main = "Posterior predictive density",
       xlab = paste("Test[", i, "]", sep = ""),  ylab = "Density",)
  plot(density(PPD.Test[, i + 1]), main = "Posterior predictive density",
       xlab = paste("Test[", i + 1, "]", sep = ""), ylab = "Density")}


M1.Y.pred <- M1.beta0.mean + M1.beta1.mean * X.test
M2.Y.pred <- M2.beta0.mean + M2.beta1.mean * X.test^3
Y.pred <- M1.weight * M1.Y.pred + M2.weight * M2.Y.pred

cat(sqrt(mean((Y.pred - Y.test)^2)))

par(mfrow=c(1,1))
ts.plot(Y.pred, col='red', ylab="Y", xlab="index",
        ylim=c(min(min(Y.pred), min(Y.test)), max(max(Y.pred), max(Y.test))))
lines(Y.test, col='blue')
legend("topleft", legend = c("Y.pred", "Y.test"),
       col = c("red", "blue"), lty = 1, cex = 1)

