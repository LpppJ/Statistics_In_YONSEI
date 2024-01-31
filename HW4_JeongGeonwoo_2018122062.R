# STA3127 Statistical Computing
# HW4
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
dev.off()
set.seed(2018122062)

##### Q1 #####
### i) ###

AP <- function(x, alpha, beta) {
  delta <- 2 * alpha^beta * (1 - alpha)^beta / (alpha^beta + (1 - alpha)^beta)
  
  ifelse(x <= 0, 
         (delta^(1/beta)/gamma(1 + 1/beta)) * exp(-delta/alpha^beta * abs(x)^beta), 
         (delta^(1/beta)/gamma(1 + 1/beta)) * exp(-delta/(1 - alpha)^beta * abs(x)^beta)
  )
}

x_range <- seq(-5, 5, length.out = 1000)

params <- list(c(0.2, 0.5), c(0.8, 0.5), c(0.2, 1.5), c(0.8, 1.5))
colors <- c("red", "blue", "green", "purple")

plot(NULL, xlim=c(-5,5), ylim=c(0, 0.5), xlab="x", ylab="Density",
     main="Density of AP(alpha, beta)")

for (i in 1:length(params)) {
  lines(x_range, sapply(x_range, AP, alpha=params[[i]][1], beta=params[[i]][2]),
        col=colors[i], lwd=2)
}

legend_labels <- sprintf("alpha=%.1f, beta=%.1f",
                         sapply(params, `[[`, 1), sapply(params, `[[`, 2))
legend("topright", legend=legend_labels, fill=colors)

### iii) ###
APAP <- function(x, alpha, beta, gamma){
  return (AP(x,alpha,beta) / AP(x,gamma,1))
}

params <- list(c(0.5, 0.5, 0.3), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.7))
colors <- c("red", "blue", "green")
plot(NULL, xlim=c(-5,5), ylim=c(0, 3), xlab="x", ylab="Density",
     main="Density of AP(alpha, beta)/AP(gamma, 1)")
for (i in 1:length(params)) {
  lines(x_range, sapply(x_range, APAP, alpha=params[[i]][1], beta=params[[i]][2],
                        gamma=params[[i]][3]), col=colors[i], lwd=2)
}
legend_labels <- sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f",
                         sapply(params, `[[`, 1), sapply(params, `[[`, 2), sapply(params, `[[`, 3))
legend("topright", legend=legend_labels, fill=colors)

### iv) ###
params <- list(c(0.5, 1.5, 0.3), c(0.5, 1.5, 0.5), c(0.5, 1.5, 0.7))
colors <- c("red", "blue", "green")
plot(NULL, xlim=c(-5,5), ylim=c(0, 3), xlab="x", ylab="Density",
     main="Density of AP(alpha, beta)/AP(gamma, 1)")
for (i in 1:length(params)) {
  lines(x_range, sapply(x_range, APAP, alpha=params[[i]][1], beta=params[[i]][2],
                        gamma=params[[i]][3]), col=colors[i], lwd=2)
}
legend_labels <- sprintf("alpha=%.1f, beta=%.1f, gamma=%.1f",
                         sapply(params, `[[`, 1), sapply(params, `[[`, 2), sapply(params, `[[`, 3))
legend("topright", legend=legend_labels, fill=colors)

### v) ###
alpha = 0.8
beta = 1.5

best_c <- function(alpha, beta){
  delta <- 2 * alpha^beta * (1 - alpha)^beta / (alpha^beta + (1 - alpha)^beta)
  g_constant <- gamma(1 + 1/beta) * 2*alpha*(1-alpha)
  exp_part1 <- (2/beta * alpha - 2*alpha)
  exp_part2 <- alpha * (alpha^beta + (1-alpha)^beta) / alpha^beta / beta
  return (delta^(1/beta) / g_constant * exp(-exp_part1*exp_part2^(1/(beta-1))))
}

envelope <- function(x, alpha, beta){
  return (best_c(alpha, beta) * AP(x, alpha, beta=1))
}

colors <- c("red", "blue")
plot(NULL, xlim=c(-5,5), ylim=c(0, 0.5), xlab="x", ylab="Density",
     main="Density of AP(alpha, beta) and cAP(alpha, 1)")
lines(x_range, sapply(x_range, AP, alpha=alpha, beta=beta), col=colors[1], lwd=2)
lines(x_range, sapply(x_range, envelope, alpha=alpha, beta=beta), col=colors[2], lwd=2)
legend_labels <- c("Target : AP(0.8, 1.5)", "Envelope : cAP(0.8, 1)")
legend("topright", legend=legend_labels, fill=colors)

### vi) ###
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

n_samples = 10^5
samples_target = rep(NA, n_samples)
now_index = 1
reject = 0
bestc = best_c(alpha, beta)

while (is.na(samples_target[n_samples])){
  sample_proposal <- AL_cdf_inv(runif(1), m=0, lambda=2*sqrt(alpha*(1-alpha)),
                                kappa=sqrt(alpha/(1-alpha)))
  if (runif(1) <= AP(sample_proposal, alpha=alpha, beta=beta)
      / (bestc * AP(sample_proposal, alpha=alpha, beta=1))){
    samples_target[now_index] = sample_proposal
    now_index = now_index + 1
  } else {
    reject = reject + 1
  }
}

hist(samples_target, xlim=c(-10,5), ylim=c(0,0.4), breaks = 100, freq=FALSE,
     main="Histogram of samples of target AP(0.8, 1.5)", xlab="sample")
lines(seq(-10, 5, length.out = 10000), sapply(seq(-10, 5, length.out = 10000),
                                              AP, alpha=alpha, beta=beta), col=colors[1] ,lwd=2)
lines(seq(-10, 5, length.out = 10000), sapply(seq(-10, 5, length.out = 10000),
                                              envelope, alpha=alpha, beta=beta), col=colors[2] ,lwd=2)
legend_labels <- c("Target : AP(0.8, 1.5)", "Envelope : cAP(0.8, 1)")
legend("topleft", legend=legend_labels, fill=colors)

cat("Acceptance probability of the algorithm:", n_samples / (n_samples+reject),"\n",
    "True acceptance probability:", 1/best_c(0.8, 1.5))

##### Q2 #####

rm(list=ls())
dev.off()
set.seed(2018122062)

### iii) ###
x_range <- seq(-5, 5, length.out = 1000)

g <- function(x, mu, sigma){
  # dnorm
  return (1/sqrt(2*pi)/sigma * exp(-(x-mu)^2 / 2*sigma^2))
}

f <- function(x, a, b, mu, sigma){
  # truncated dnorm
  dnorm_part <- (1/sqrt(2*pi)/sigma * exp(-(x-mu)^2 / 2*sigma^2))
  normalizing <- 1 - pnorm(b,mu,sigma) + pnorm(a,mu,sigma)
  indicator <- 1*(x < a) + 1*(x > b)
  return (dnorm_part / normalizing * indicator)
}

envelop <- function(x, mu, sigma){
  return (c_ * g(x,mu,sigma))
}

n_samples = 10^5
samples_target = rep(NA, n_samples)
now_index = 1
reject = 0
c_ = 1 / (1 - pnorm(2, 1, 1) + pnorm(0.5, 1, 1))

while (is.na(samples_target[n_samples])){
  sample_proposal <- qnorm(runif(1), 1, 1)
  if (runif(1) <= f(sample_proposal, 0.5, 2, 1, 1)
      / g(sample_proposal, 1, 1) / c_){
    samples_target[now_index] = sample_proposal
    now_index = now_index + 1
  } else {
    reject = reject + 1
  }
}

cat("Estimated acceptance probability:",
    length(samples_target)/(length(samples_target)+reject),"\n",
    "True acceptance probability:",1/c_)

colors <-("red")
hist(samples_target, ylim=c(0,1), breaks = 100, freq = FALSE,
     main="Histogram of samples of target TN", xlab="sample")
lines(x_range, sapply(x_range, f, a=0.5, b=2, mu=1, sigma=1),
      col=colors[1], lwd=2)
legend_labels <- c("Target : TN")
legend("topright", legend=legend_labels, fill=colors)

### iv) ###
x_range <- seq(-5, 7, length.out = 10000)

envelop <- function(x, mu, sigma){
  return (c_ * g(x,mu,sigma))
}

n_samples = 10 # Even sampling 10 takes a few tens of seconds.
samples_target = rep(NA, n_samples)
now_index = 1
reject = 0
c_ = 1 / (1 - pnorm(5.5, 1, 1) + pnorm(-4, 1, 1))

while (is.na(samples_target[n_samples])){
  sample_proposal <- qnorm(runif(1), 1, 1)
  if (runif(1) <= f(sample_proposal, -4, 5.5, 1, 1)
      / g(sample_proposal, 1, 1) / c_){
    samples_target[now_index] = sample_proposal
    now_index = now_index + 1
  } else {
    reject = reject + 1
  }
}

colors <- c("red", "blue")
hist(samples_target, xlim = c(-5,7), ylim=c(0,10), breaks = 100,
     freq = FALSE, main="Histogram of samples of target TN", xlab="sample")
lines(x_range, sapply(x_range, f, a=-4, b=5.5, mu=1, sigma=1),
      col=colors[1], lwd=2)
lines(x_range, sapply(x_range, envelop, mu=1, sigma=1),
      col=colors[2], lwd=2)
legend_labels <- c("Target : TN", "Envelope : cN(1,1)")
legend("top", legend=legend_labels, fill=colors)


### vi) ###

x_range <- seq(-7, 7, length.out = 10000)

g <- function(U, a, b, mu, sigma) {
  pnorm1 <- pnorm(a, mu, sigma)
  pnorm2 <- pnorm(b, mu, sigma)
  pnorm12 <- 1-pnorm2+pnorm1
  if (U <= pnorm1 / pnorm12) {
    return(qnorm(U * pnorm12, 1, 1))
  } else {
    return(qnorm(pnorm12*U +pnorm2 - pnorm1, 1, 1))
  }
}
n_samples = 10^5
samples_target = rep(NA, n_samples)
U <- runif(n_samples)
for (i in 1:n_samples){
  samples_target[i] <- g(U[i], a=-4, b=5.5, mu=1, sigma=1)
}

colors <- c("red")
hist(samples_target, xlim = c(-7,7), ylim=c(0,5), breaks = 100,
     freq = FALSE, main="Histogram of samples of target TN", xlab="sample")
lines(x_range, sapply(x_range, f, a=-4, b=5.5, mu=1, sigma=1),
      col=colors[1], lwd=2)
legend_labels <- c("Target : TN")
legend("top", legend=legend_labels, fill=colors)

