# STA3105 Bayesian Statistics
# HW5
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)

library(mvtnorm);library(fields);library(classInt);library(MCMCpack)
###### Q1 ######
### (a) ###
n = 300000
X <- mvrnorm(n =n, mu = rep(0, 2), Sigma = diag(2))

### (b) ###
beta.true <- c(0.5, 1, 2, -1)
X_design <- cbind( rep(1,n), X, X[,1]*X[,2] )
Y <- X_design %*% beta.true + rnorm(n, 0, 1)

##### Q2 #####
p = ncol(X_design)
m.beta <- 0; v.beta <- 10
a.s2 <- 0.01; b.s2 <- 0.01

B = 1000

# sample space (with initialization)
beta.samps <- matrix(NA, nrow = p, ncol = B)
beta.samps[,1] <- rep(1,p)
s2.samps <- matrix(NA, nrow = 1, ncol = B)
s2.samps[1] <- 1

# Gibbs sampler
ptm <- proc.time()
for(i in 2:B){
  
  ## beta[i] | s2[i-1]
  V <- solve( t(X_design) %*% X_design / s2.samps[i-1] + diag(p) / v.beta )
  m <- V %*% ( t(X_design) %*% Y / s2.samps[i-1] )
  beta.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  # s2[i] | beta[i]
  a <- n/2 + a.s2
  b <- b.s2 + t( Y - X_design %*% beta.samps[,i] ) %*% ( Y - X_design %*% beta.samps[,i] ) / 2
  s2.samps[i] <- rinvgamma(1, shape=a, scale=b)
  s2.samps <- t(s2.samps)
}
Rtime <- proc.time() - ptm

# ts plot
ts.plot(beta.samps[1,], main="Density Plot of beta0", xlab="iter", ylab="beta")
ts.plot(beta.samps[2,], main="Density Plot of beta1", xlab="iter", ylab="beta")
ts.plot(beta.samps[3,], main="Density Plot of beta2", xlab="iter", ylab="beta")
ts.plot(beta.samps[4,], main="Density Plot of beta3", xlab="iter", ylab="beta")
ts.plot(s2.samps, main="Density Plot of sigma2", xlab="iter", ylab="sigma2")

# burn-in
beta.samps <- beta.samps[,51:B]
s2.samps <- s2.samps[51:B]


# density plot
plot(density(beta.samps[1,]), main="Density Plot of beta0", xlab="beta", ylab="Density", lwd=2)
plot(density(beta.samps[2,]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2)
plot(density(beta.samps[3,]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2)
plot(density(beta.samps[4,]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2)
plot(density(s2.samps), main="Density Plot of sigma2", xlab="sigma2", ylab="Density", lwd=2)

# 95% HPD intervals
cat("95% HPD intervals of beta0:", HPDinterval(as.mcmc(beta.samps[1,]), prob = 0.95), "\n",
    "95% HPD intervals of beta1:", HPDinterval(as.mcmc(beta.samps[2,]), prob = 0.95), "\n",
    "95% HPD intervals of beta2:", HPDinterval(as.mcmc(beta.samps[3,]), prob = 0.95), "\n",
    "95% HPD intervals of beta3:", HPDinterval(as.mcmc(beta.samps[4,]), prob = 0.95), "\n",
    "95% HPD intervals of sigma2:", HPDinterval(as.mcmc(s2.samps), prob = 0.95) )

# posterior mean
cat("posterior mean of beta0:", mean(beta.samps[1,]), "\n",
    "posterior mean of beta1:", mean(beta.samps[2,]), "\n",
    "posterior mean of beta2:", mean(beta.samps[3,]), "\n",
    "posterior mean of beta3:", mean(beta.samps[4,]), "\n",
    "posterior mean of sigma2:", mean(s2.samps) )
    

# acceptance probability
# (Because I used Gibbs sampler, samples was always accepted !)
cat("Acceptance probability of beta0 :", length(unique(beta.samps[1,]))/B ,"\n",
    "Acceptance probability of beta1 :", length(unique(beta.samps[2,]))/B ,"\n",
    "Acceptance probability of beta2 :", length(unique(beta.samps[3,]))/B ,"\n",
    "Acceptance probability of beta3 :", length(unique(beta.samps[4,]))/B ,"\n",
    "Acceptance probability of sigma2 :", length(unique(s2.samps))/B )

# effective sample size
cat("Effective Sample size of beta0:", effectiveSize(beta.samps[1,]), "\n",
    "Effective Sample size of beta1:", effectiveSize(beta.samps[2,]), "\n",
    "Effective Sample size of beta2:", effectiveSize(beta.samps[3,]), "\n",
    "Effective Sample size of beta3:", effectiveSize(beta.samps[4,]), "\n",
    "Effective Sample size of sigma2:", effectiveSize(s2.samps) )
    
##### Q3 #####
library(Rcpp); library(RcppArmadillo)
setwd("/Users/gunwoojung/Desktop/R_wd/Bayesian_Stat")
sourceCpp("HW5_JeongGeonwoo_2018122062.cpp")   

ptm <- proc.time()
MCMC_Rcpp <- RcppGibbs(B, X_design, Y,  m.beta, v.beta, a.s2, b.s2)
Rcpptime <- proc.time() - ptm
beta_samps <- MCMC_Rcpp$beta_samps
s2_samps <- MCMC_Rcpp$s2_samps

# ts plot
ts.plot(beta_samps[1,], main="Density Plot of beta0", xlab="iter", ylab="beta")
ts.plot(beta_samps[2,], main="Density Plot of beta1", xlab="iter", ylab="beta")
ts.plot(beta_samps[3,], main="Density Plot of beta2", xlab="iter", ylab="beta")
ts.plot(beta_samps[4,], main="Density Plot of beta3", xlab="iter", ylab="beta")
ts.plot(s2_samps, main="Density Plot of sigma2", xlab="iter", ylab="sigma2")

# burn-in
beta_samps <- beta_samps[,51:B]
s2_samps <- s2_samps[51:B]

# density plot
plot(density(beta_samps[1,]), main="Density Plot of beta0", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samps[2,]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samps[3,]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2)
plot(density(beta_samps[4,]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2)
plot(density(s2_samps), main="Density Plot of sigma2", xlab="sigma2", ylab="Density", lwd=2)

# 95% HPD intervals
cat("95% HPD intervals of beta0:", HPDinterval(as.mcmc(beta_samps[1,]), prob = 0.95), "\n",
    "95% HPD intervals of beta1:", HPDinterval(as.mcmc(beta_samps[2,]), prob = 0.95), "\n",
    "95% HPD intervals of beta2:", HPDinterval(as.mcmc(beta_samps[3,]), prob = 0.95), "\n",
    "95% HPD intervals of beta3:", HPDinterval(as.mcmc(beta_samps[4,]), prob = 0.95), "\n",
    "95% HPD intervals of sigma2:", HPDinterval(as.mcmc(s2_samps), prob = 0.95) )

# posterior mean
cat("posterior mean of beta0:", mean(beta_samps[1,]), "\n",
    "posterior mean of beta1:", mean(beta_samps[2,]), "\n",
    "posterior mean of beta2:", mean(beta_samps[3,]), "\n",
    "posterior mean of beta3:", mean(beta_samps[4,]), "\n",
    "posterior mean of sigma2:", mean(s2_samps) )


# acceptance probability
# (Because I used Gibbs sampler, samples was always accepted !)
cat("Acceptance probability of beta0 :", length(unique(beta_samps[1,]))/B ,"\n",
    "Acceptance probability of beta1 :", length(unique(beta_samps[2,]))/B ,"\n",
    "Acceptance probability of beta2 :", length(unique(beta_samps[3,]))/B ,"\n",
    "Acceptance probability of beta3 :", length(unique(beta_samps[4,]))/B ,"\n",
    "Acceptance probability of sigma2 :", length(unique(s2_samps))/B )

# effective sample size
cat("Effective Sample size of beta0:", effectiveSize(beta_samps[1,]), "\n",
    "Effective Sample size of beta1:", effectiveSize(beta_samps[2,]), "\n",
    "Effective Sample size of beta2:", effectiveSize(beta_samps[3,]), "\n",
    "Effective Sample size of beta3:", effectiveSize(beta_samps[4,]), "\n",
    "Effective Sample size of sigma2:", effectiveSize(s2_samps) )

plot(density(beta.samps[1,]), main="Density Plot of beta0", xlab="beta", ylab="Density", lwd=2, col="blue")
lines(density(beta_samps[1,]), main="Density Plot of beta0", xlab="beta", ylab="Density", lwd=2, col="red")
legend("topright", legend=c("R","Rcpp"), fill=c("blue","red"))

plot(density(beta.samps[2,]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2, col="blue")
lines(density(beta_samps[2,]), main="Density Plot of beta1", xlab="beta", ylab="Density", lwd=2, col="red")
legend("topright", legend=c("R","Rcpp"), fill=c("blue","red"))

plot(density(beta.samps[3,]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2, col="blue")
lines(density(beta_samps[3,]), main="Density Plot of beta2", xlab="beta", ylab="Density", lwd=2, col="red")
legend("topright", legend=c("R","Rcpp"), fill=c("blue","red"))

plot(density(beta.samps[4,]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2, col="blue")
lines(density(beta_samps[4,]), main="Density Plot of beta3", xlab="beta", ylab="Density", lwd=2, col="red")
legend("topright", legend=c("R","Rcpp"), fill=c("blue","red"))

plot(density(s2.samps), main="Density Plot of sigma2", xlab="sigma2", ylab="Density", lwd=2, col="blue")
lines(density(s2_samps), main="Density Plot of sigma2", xlab="sigma2", ylab="Density", lwd=2, col="red")
legend("topright", legend=c("R","Rcpp"), fill=c("blue","red"))

##### Q4 #####
library(parallel) # for checking number of cores (not for implementation)
numcore <- detectCores(); numcore

# this is for openmp 
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

ptm <- proc.time()
MCMC_RcppParallel <- RcppGibbsparallel(B, X_design, Y,  m.beta, v.beta, a.s2, b.s2, numcore)
RcppParalleltime <- proc.time() - ptm

beta_samps_ <- MCMC_RcppParallel$beta_samps
s2_samps_ <- matrix(MCMC_RcppParallel$s2_samps, ncol=B, nrow=numcore)

Rtime; Rcpptime; RcppParalleltime

# burn in
burn_in <- 50
beta_samps_ <- beta_samps_[,(burn_in+1):B,]
s2_samps_ <- s2_samps_[,(burn_in+1):B]

# ts plot
ts.plot(beta_samps_[,,1][1,],main="Trace Plot of beta0", xlab="iter", ylab="beta0")
for(i in 2:numcore){ lines(beta_samps_[,,i][1,],col=i) }
ts.plot(beta_samps_[,,1][2,],main="Trace Plot of beta1", xlab="iter", ylab="beta1")
for(i in 2:numcore){ lines(beta_samps_[,,i][2,],col=i) }
ts.plot(beta_samps_[,,1][2,],main="Trace Plot of beta2", xlab="iter", ylab="beta2")
for(i in 2:numcore){ lines(beta_samps_[,,i][2,],col=i) }
ts.plot(beta_samps_[,,1][2,],main="Trace Plot of beta3", xlab="iter", ylab="beta3")
for(i in 2:numcore){ lines(beta_samps_[,,i][2,],col=i) }
ts.plot(s2_samps_[1,],main="Trace Plot of sigma2", xlab="iter", ylab="sigma2")
for(i in 2:numcore){ lines(s2_samps_[i,],col=i) }

# density plot
plot(density(beta_samps_[,,1][1,]), main="Density Plot of beta0", xlab="beta0", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][1,]), col=i) }
plot(density(beta_samps_[,,1][2,]), main="Density Plot of beta1", xlab="beta1", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][2,]), col=i) }
plot(density(beta_samps_[,,1][3,]), main="Density Plot of beta2", xlab="beta2", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][3,]), col=i) }
plot(density(beta_samps_[,,1][4,]), main="Density Plot of beta3", xlab="beta3", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][4,]), col=i) }
plot(density(s2_samps_[1,]), main="Density Plot of sigma2", xlab="sigma2", ylab="Density")
for(i in 2:numcore){ lines(density(s2_samps_[i,]), col=i) }

# 95% HPD intervals
print_interval <- function(core){
  for (j in 1:p){
    for (i in 1:core){
      cat("95% HPD intervals of beta",j-1,"/shard", i,":",
          HPDinterval(as.mcmc(beta_samps_[,,i][j,]), prob = 0.95), "\n")}
    cat("\n")}
  for (i in 1: core){
    cat("95% HPD intervals of sigma2 /shard", i,":",
        HPDinterval(as.mcmc(s2_samps_[i,]), prob = 0.95), "\n")}}
print_interval(numcore)


# posterior mean
print_posterior_mean <- function(core){
  for (j in 1:p){
    for (i in 1:core){
      cat("posterior mean of beta",j-1, "/shard:", i,":",
          mean(beta_samps_[,,i][j,]), "\n")}
    cat("\n")}
  for (i in 1: core){
    cat("posterior mean of sigma2 /shard", i,":",
        mean(s2_samps_[i,], prob = 0.95), "\n")}}

print_posterior_mean(numcore)


# acceptance probability
# (Because I used Gibbs sampler, samples was always accepted !)
print_acc_prob <- function(core){
  for (j in 1:p){
    for (i in 1:core){
      cat("Acceptance probability of beta", j-1, "/shard", i,":",
          length(unique(beta_samps_[,,i][j,]))/B ,"\n")}
    cat("\n")}
  for (i in 1: core){
    cat("Acceptance probability of sigma2 /shard", i,":",
        length(s2_samps_[i,])/B, "\n")}}

print_acc_prob(numcore)

# effective sample size
print_ess <- function(core){
  for (j in 1:p){
    for (i in 1:core){
      cat("Effective Sample size of beta", j-1, "/shard", i,":",
          effectiveSize(beta_samps_[,,i][p,]), "\n")}
    cat("\n")}
  for (i in 1: core){
    cat("Acceptance probability of sigma2 /shard", i,":",
        effectiveSize(s2_samps_[i,]), "\n")}}

print_ess(numcore)



# (c)
# density plot
plot(density(beta_samps_[,,1][1,]), main="Density Plot of beta0", xlab="beta0", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][1,]), col=i) }
lines(density(beta_samps[1,]))

plot(density(beta_samps_[,,1][2,]), main="Density Plot of beta1", xlab="beta1", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][2,]), col=i) }
lines(density(beta_samps[2,]))

plot(density(beta_samps_[,,1][3,]), main="Density Plot of beta2", xlab="beta2", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][3,]), col=i) }
lines(density(beta_samps[3,]))

plot(density(beta_samps_[,,1][4,]), main="Density Plot of beta3", xlab="beta3", ylab="Density")
for(i in 2:numcore){ lines(density(beta_samps_[,,i][4,]), col=i) }
lines(density(beta_samps[4,]))

plot(density(s2_samps_[1,]), main="Density Plot of sigma2", xlab="sigma2", ylab="Density")
for(i in 2:numcore){ lines(density(s2_samps_[i,]), col=i) }
lines(density(s2_samps))
