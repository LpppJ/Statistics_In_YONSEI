# STA3105 Bayesian Statistics
# HW4
# Name : Jeong Geonwoo
# ID : 2018122062

set.seed(2018122062)
# rm(list=ls())
library(mvtnorm);library(fields);library(classInt);library(MCMCpack)

##### Problem 1 #####
# Generating Spatial locations
n <- 1000
spatial_locations <- matrix(runif(2 * n), ncol = 2)
X <- matrix(runif(2*n, -1, 1), ncol = 2)

# Prior layer (using true parameter values)
beta <- c(1,1)
rho <- 0.2
sigma2 <- 1

# Process layer
distances <- as.matrix(dist(spatial_locations))
Gamma <- exp(-distances/rho)
eta <- rmvnorm(1, X%*%beta, sigma2*Gamma)

# Data layer
tau2 <- 0.2
Y <- rmvnorm(1, eta, tau2*diag(n))

# plotRF
# input of plotRF function 
# (1) dat: variable to be plotted
# (2) range: variable to be plotted
# (3) label: title of the plot
# (4) location: spatial location of the variable to be plotted
plotRF <- function(dat,rangeDat,label,location,length.out=10){
  breaks <- seq(range(rangeDat,na.rm = TRUE)[1],
                range(rangeDat,na.rm = TRUE)[2],
                length.out=length.out)
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(dat, n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=location[,1],y=location[,2],col=col, pch=16,
       main=label)
}
plotRF(dat=Y, rangeDat=Y, label="2-dim plot of simulated Y", location=spatial_locations)

##### Problem 3 #####
# Prior parameters
m.beta <- rep(0, 2); V.beta <- 10 * diag(2)
a.s2 <- 0.001; b.s2 <- 0.001
a.t2 <- 0.001; b.t2 <- 0.001

rhoseq <- seq(0.001, 1, length = 100)
a.rho <- 0.01; b.rho <- 0.01
plot(rhoseq, dgamma(rhoseq, shape = a.rho, scale = b.rho)) # prior for rho

B <- 2000

beta.samps <- matrix(NA, nrow = 2, ncol = B)
beta.samps[,1] <- c(1,1)

s2.samps <- t2.samps <- rho.samps <- rep(NA, B)
s2.samps[1] <- 1
rho.samps[1] <- 0.2
t2.samps[1] <- 0.2

eta.obs.samps <- matrix(NA, nrow = n, ncol = B)
eta.obs.samps[,1] = X%*%beta.samps[,1]

# the variance of proposed rho in following MCMC
v.prop <- 0.01

# Defien Gamma
Gamma <- exp(-distances/rho.samps[1])
Ginv <- solve(Gamma)
# For vectorized computing, I converted Y to matrix
Y <- matrix(Y, ncol=1)

# We will use the first 500 obs in Problem 3 !
# "name_" means that the first 500 value of "name"
X_ <- X[1:(n/2),] # (note!) X[1:(n/2),] and X[1:n/2,] is different....
Y_ <- Y[1:(n/2),]
distances_ <- distances[1:(n/2),1:(n/2)]
Gamma_ <- exp( -distances_ / rho.samps[1] )
Ginv_ <- solve(Gamma_)
eta.obs.samps_ <- eta.obs.samps[1:(n/2),]
##

i=2
# Gibbs sampler (MH alghrithm for rho)
for(i in 2:B){
  
  if(i%%100==0) print(i)
  
  ## eta_obs | Rest
  V <- solve( diag(n/2) / t2.samps[i-1] + Ginv_ / s2.samps[i-1] )
    m <- V %*% ( Y_ / t2.samps[i-1] + Ginv_ %*% X_ %*% beta.samps[,i-1] / s2.samps[i-1] )
    eta.obs.samps_[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  ## beta | Rest
  V <-  solve( diag(2) / 10 + t(X_) %*% Ginv_ %*% X_ / s2.samps[i-1])
    m <- V %*% ( t(X_) %*% Ginv_ %*% eta.obs.samps_[,i] / s2.samps[i-1])
    beta.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")

  ## s2 | Rest
  a <- a.s2 + n/4 # n/4 = (the number of obs)/2 = (n/2)/2 = 250
    b <- b.s2 + t( eta.obs.samps_[,i] - X_ %*% beta.samps[,i] ) %*% Ginv_ %*% ( eta.obs.samps_[,i] - X_ %*% beta.samps[,i] ) / 2
    s2.samps[i] <- rinvgamma(1, shape=a, scale=b)
  
  ## t2 | Rest
  a <- a.t2 + n/4 # n/4 = (the number of obs)/2 = (n/2)/2 = 250
    b <- b.t2 + t( Y_ - eta.obs.samps_[,i] ) %*% ( Y_ - eta.obs.samps_[,i] ) / 2
    t2.samps[i] <- rinvgamma(1, shape=a, scale=b)
  
  ## rho | Rest   
  # MH update
  rho_proposal <- rnorm(1, mean=rho.samps[i-1], sd=sqrt(v.prop))
  # rho should be greater than 0
  while (rho_proposal < 0){
    rho_proposal = rnorm(1, mean=rho.samps[i-1], sd=sqrt(v.prop))
  }
  
  # the function for kernel
  # input : (rho_proposal) or (rho_previous)
  # output : the value which is proportional to log_posterior
  log_rho_kernel <- function(rho){
    # Gamma (using given rho)
    G_ <- exp( -distances_ / rho ) 
    # I used (dmvnorm) instead (kernel Multivariate normal distribution PDF)
    # reason : det(Gamma) = 0
    normal_ <- dmvnorm( eta.obs.samps_[,i], mean=X_ %*% beta.samps[,i],  sigma = G_ * s2.samps[i], log = TRUE)
    # I used kernel of Gamma distribution PDF
    gamma_ <- (a.rho - 1) * log( rho ) - b.rho * rho
    # gamma_ <- dgamma(rho, shape=a.rho, scale=b.rho)
    # log(AB) = logA + logB
    return (normal_ + gamma_)
  }
  # log(A/B) = logA - logB
  acc_prob <- log_rho_kernel(rho_proposal) - log_rho_kernel(rho.samps[i-1])
  if ( acc_prob > log(runif(1)) ){
    rho.samps[i] <- rho_proposal
  } else {
    rho.samps[i] <- rho.samps[i-1]
  }
  # Update Gamma and Ginv (using new rho)
  Gamma_ <- exp(-distances_ / rho.samps[i])
  Ginv_ <- solve(Gamma_)
}

# Acceptance rate of rho
length(unique(rho.samps))/B

# effective sample size
cat("Effective Sample size of beta.samps[1]:", effectiveSize(beta.samps[1,]))
cat("Effective Sample size of beta.samps[2]:", effectiveSize(beta.samps[2,]))
cat("Effective Sample size of eta.obs.samps:", effectiveSize(eta.obs.samps_[,1]))
cat("Effective Sample size of s2.samps:", effectiveSize(s2.samps))
cat("Effective Sample size of t2.samps:", effectiveSize(t2.samps))
cat("Effective Sample size of rho.samps:", effectiveSize(rho.samps))

# trace plot (Need burn-in ?)
# par(mfrow=c(3,2))
ts.plot(beta.samps[1,])
ts.plot(beta.samps[2,])
ts.plot(eta.obs.samps_[1,])
ts.plot(s2.samps)
ts.plot(t2.samps)
ts.plot(rho.samps)

# ACF plot (Need thinning ?)
# par(mfrow=c(3,2))
acf(beta.samps[1,])
acf(beta.samps[2,])
acf(eta.obs.samps_[,1])
acf(s2.samps)
acf(t2.samps)
acf(rho.samps)

# posterior mean
cat("Posterior mean of beta.samps[1,] :", mean(beta.samps[1,]))
cat("Posterior mean of beta.samps[2,] :", mean(beta.samps[2,]))
cat("Posterior mean of eta.obs.samps_[1,] :", mean(eta.obs.samps_[1,]))
cat("Posterior mean of s2.samps :", mean(s2.samps))
cat("Posterior mean of t2.samps :", mean(t2.samps))
cat("Posterior mean of rho.samps :", mean(rho.samps))

# 95% HPD Interval
cat("95% HPD interval of beta.samps[1,]:", HPDinterval(as.mcmc(beta.samps[1,])))
cat("95% HPD interval of beta.samps[2,]:", HPDinterval(as.mcmc(beta.samps[2,])))
cat("95% HPD interval of eta.obs.samps_[1,]:", HPDinterval(as.mcmc(eta.obs.samps_[1,])))
cat("95% HPD interval of s2.samps:", HPDinterval(as.mcmc(s2.samps)))
cat("95% HPD interval of t2.samps:", HPDinterval(as.mcmc(t2.samps)))
cat("95% HPD interval of rho.samps:", HPDinterval(as.mcmc(rho.samps)))

# After MCMC Processing
# Burn-in and Thinning is Needed !

##### Problem 4 #####
# Prediction

X.pred <- X[(n/2+1):n,]

dcross <- rdist(spatial_locations[1:(n/2),], spatial_locations[(n/2+1):n,])
dpred <- rdist(spatial_locations[(n/2+1):n,])

eta.pred <- matrix(NA, nrow = nrow(X.pred), ncol = B)

j=1

for(j in 1:B){
  
  if(j%%100==0) print(j)
  
  # update gamma
  Gamma_ = exp( -distances_/rho.samps[j] )
  Ginv_ = solve(Gamma_)
  Gamma.pred <- exp( -dpred / rho.samps[j] )
  Gamma.cross <- exp( -dcross / rho.samps[j] )
  
  m <- X.pred %*% beta.samps[,j]
        + t(Gamma.cross) %*% Ginv_ %*% (matrix(eta.obs.samps_[,i], ncol=1) - X_%*% beta.samps[,j])
  V <- s2.samps[j] * (Gamma.pred - t(Gamma.cross) %*% Ginv_ %*% Gamma.cross)
  eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}

eta.predicted <- rbind(eta.obs.samps_, eta.pred)
eta.pred_post = rep(NA, dim(eta.pred)[1])
eta.pred_sd = rep(NA, dim(eta.pred)[1])
for (i in 1:dim(eta.predicted)[1]){
  eta.pred_post[i] = mean(eta.predicted[i,])
  eta.pred_sd[i] = sd(eta.predicted[i,])
}

plotRF(dat=eta, rangeDat=eta, label="Surface with the true eta", location=spatial_locations)
plotRF(dat=eta.pred_post, rangeDat=eta.pred_post, label="Surface with the eta.pred_post", location=spatial_locations)
plotRF(dat=eta.pred_sd, rangeDat=eta.pred_sd, label="Surface with the posterior standard deviation", location=spatial_locations)

