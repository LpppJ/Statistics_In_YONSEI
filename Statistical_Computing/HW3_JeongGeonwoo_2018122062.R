# STA3127 Statistical Computing
# HW3
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)

##### Q1 #####
AL_pdf <- function(x, m, lambda, kappa) {
  pdf <- (lambda / (kappa + 1/kappa)) * exp(-lambda * (x-m) * sign(x-m) * kappa^{sign(x-m)})
  return (pdf)
}

AL_cdf_inv <- function(U, m, lambda, kappa) {
  sapply(U, function(Ui) {
    if (Ui < kappa^2 / (1 + kappa^2)) {
      cdf_inv <- kappa / lambda * log((1 + kappa^2) / kappa^2 * Ui) + m
    } else {
      cdf_inv <- - 1 / (kappa * lambda) * log((1 + kappa^2)*(1-Ui)) + m
    }
    return(cdf_inv)
  })
}

U <- runif(10^6)
samples <- AL_cdf_inv(U, m=1, lambda=1, kappa=2)
hist(samples, breaks = 100, freq=FALSE)
curve(AL_pdf(x, m=1, lambda=1, kappa=2), from = min(samples), to = max(samples), add=TRUE, lwd=2)

##### Q2 #####
Geometric_cdf <- function(x, q) {
  sapply(x, function(xi) {
    cdf <- 1-(1-q)^{max(0, floor(xi))}
    return(cdf)
  })
}

Geometric_cdf_inv <- function(x, q) {
  return(ceiling(log(1-x) / log(1-q)))
}

p = 0.7
q = 1-p
r = 10

n_samples=10^6
U = runif(n_samples,0,1)
samples <- rep(NA, n_samples)


for (i in 1:n_samples){
  sum_geometric <- sum(Geometric_cdf_inv(runif(r), q))
  neg_binorm <- sum_geometric - r
  samples_i <- sin(neg_binorm)
  samples[i] <- samples_i
}

mean(samples)