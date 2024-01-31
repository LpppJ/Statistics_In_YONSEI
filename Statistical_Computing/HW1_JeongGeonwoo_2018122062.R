# STA3127 Statistical Computing
# HW1
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())

# Parameterization
mu1 <- c(1, 0.5)
mu2 <- c(-1, -1)
sigma1 <- matrix(c(1, 0.5, 0.5, 1.5), nrow = 2)
sigma2 <- matrix(c(2, -0.2, -0.2, 1), nrow = 2)
weights <- c(0.4, 0.6)

############################
##### hw 1 - problem 1 #####
##### [-10, 10]^2 grid #####
############################

data.grid <- expand.grid(x.1 = seq(-10, 10, length.out = 101), x.2 = seq(-10, 10, length.out = 101))

multivariate_normal_pdf <- function(x, mean, cov) {
  d <- length(mean)
  det_cov <- det(cov)
  inv_cov <- solve(cov)
  constant <- 1 / ((2 * pi)^(d/2) * sqrt(det_cov))
  exponent <- -0.5 * t(x - mean) %*% inv_cov %*% (x - mean)
  pdf_value <- constant * exp(exponent)
  return (pdf_value)
}

probabilities1 <- apply(data.grid, 1, function(x) {
  multivariate_normal_pdf(x, mu1, sigma1)
})
probabilities2 <- apply(data.grid, 1, function(x) {
  multivariate_normal_pdf(x, mu2, sigma2)
})

probabilities = weights[1] * probabilities1 + weights[2] * probabilities2
log_probabilities = log(probabilities)

q.samp <- cbind(data.grid, prob = log_probabilities)
contour(x = unique(q.samp$x.1), y = unique(q.samp$x.2), z = matrix(q.samp$prob, ncol = 101),
        main = "Contour Plot of Gaussian Mixture Distribution",
        xlab = "x1", ylab = "x2")

############################
##### hw 1 - problem 1 #####
##### [-50, 50]^2 grid #####
############################

data.grid <- expand.grid(x.1 = seq(-50, 50, length.out = 101), x.2 = seq(-50, 50, length.out = 101))

multivariate_normal_pdf <- function(x, mean, cov) {
  d <- length(mean)
  det_cov <- det(cov)
  inv_cov <- solve(cov)
  constant <- 1 / ((2 * pi)^(d/2) * sqrt(det_cov))
  exponent <- -0.5 * t(x - mean) %*% inv_cov %*% (x - mean)
  pdf_value <- constant * exp(exponent)
  return (pdf_value)
}

probabilities1 <- apply(data.grid, 1, function(x) {
  multivariate_normal_pdf(x, mu1, sigma1)
})
probabilities2 <- apply(data.grid, 1, function(x) {
  multivariate_normal_pdf(x, mu2, sigma2)
})

probabilities = weights[1] * probabilities1 + weights[2] * probabilities2
log_probabilities = log(probabilities)

q.samp <- cbind(data.grid, prob = log_probabilities)
contour(x = unique(q.samp$x.1), y = unique(q.samp$x.2), z = matrix(q.samp$prob, ncol = 101),
        main = "Contour Plot of Gaussian Mixture Distribution",
        xlab = "x1", ylab = "x2")

############################
##### hw 1 - problem 1 #####
##### [-10, 10]^2 grid #####
############################
data.grid <- expand.grid(x.1 = seq(-10, 10, length.out = 101), x.2 = seq(-10, 10, length.out = 101))

log_multivariate_normal_pdf <- function(x, mean, conv) {
 d <- length(mean)
 det_cov <- det(conv)
 inv_cov <- solve(conv)
 constant <- -d/2 * log(2*pi) - 1/2 * log(det_cov)
 exponent <- -0.5 * t(x - mean) %*% inv_cov %*% (x - mean)
 log_pdf_value <- constant + exponent
 return (log_pdf_value)
}

log_prob1 <- apply(data.grid, 1, function(x) {
  log_multivariate_normal_pdf(x, mu1, sigma1)
})
log_prob2 <- apply(data.grid, 1, function(x) {
  log_multivariate_normal_pdf(x, mu2, sigma2)
})

log_weighted1 <- log(weights[1]) + log_prob1
log_weighted2 <- log(weights[2]) + log_prob2

log_prob <- numeric(length(log_weighted1))
for (i in 1:length(log_weighted1)) {
  if (log_weighted2[i] > log_weighted1[i]) {
    log_prob[i] <- log_weighted2[i] + log(1 + exp(log_weighted1[i] - log_weighted2[i]))
  } else {
    log_prob[i] <- log_weighted1[i] + log(1 + exp(log_weighted2[i] - log_weighted1[i]))
  }
}

q.samp <- cbind(data.grid, prob = log_prob)
contour(x = unique(q.samp$x.1), y = unique(q.samp$x.2), z = matrix(q.samp$prob, ncol = 101),
        main = "Contour Plot of Gaussian Mixture Distribution",
        xlab = "x1", ylab = "x2")

############################
##### hw 1 - problem 1 #####
##### [-50, 50]^2 grid #####
############################
data.grid <- expand.grid(x.1 = seq(-50, 50, length.out = 101), x.2 = seq(-50, 50, length.out = 101))

log_multivariate_normal_pdf <- function(x, mean, conv) {
  d <- length(mean)
  det_cov <- det(conv)
  inv_cov <- solve(conv)
  constant <- -d/2 * log(2*pi) - 1/2 * log(det_cov)
  exponent <- -0.5 * t(x - mean) %*% inv_cov %*% (x - mean)
  log_pdf_value <- constant + exponent
  return (log_pdf_value)
}

log_prob1 <- apply(data.grid, 1, function(x) {
  log_multivariate_normal_pdf(x, mu1, sigma1)
})
log_prob2 <- apply(data.grid, 1, function(x) {
  log_multivariate_normal_pdf(x, mu2, sigma2)
})

log_weighted1 <- log(weights[1]) + log_prob1
log_weighted2 <- log(weights[2]) + log_prob2

log_prob <- numeric(length(log_weighted1))
for (i in 1:length(log_weighted1)) {
  if (log_weighted2[i] > log_weighted1[i]) {
    log_prob[i] <- log_weighted2[i] + log(1 + exp(log_weighted1[i] - log_weighted2[i]))
  } else {
    log_prob[i] <- log_weighted1[i] + log(1 + exp(log_weighted2[i] - log_weighted1[i]))
  }
}

q.samp <- cbind(data.grid, prob = log_prob)
contour(x = unique(q.samp$x.1), y = unique(q.samp$x.2), z = matrix(q.samp$prob, ncol = 101),
        main = "Contour Plot of Gaussian Mixture Distribution",
        xlab = "x1", ylab = "x2")


############################
##### hw 1 - problem 2 #####
############################
set.seed(2018122062)
id = 18122062

# i)
options(digits=15)
scaling_factor <- 1e08

id_demical= id / scaling_factor
id_demical


middle8 <- function(number){
  # input : 16-digit number
  # output : 8-digit number which is middle of input
  result = (number %% 1e12 - number %% 1e04) / 1e04
  return (result)
  # ex.
  # middle8(1234567890123456)
  # = (0000567890123456 - 3456) / 1e04
  # = 567890120000 / 1e04
  # = 56789012
}

Middle_Square <- function(sample_size){
  samples = c()
  for (i in 1:sample_size) {
    if (i == 1) {
      samples <- c(samples, middle8(id^2))
    } else {
      samples <- c(samples, middle8((samples[length(samples)])^2))
    }
  }
  samples <- samples / scaling_factor
  return (samples)
}

Sample100 = Middle_Square(100)
ts.plot(Sample100, main = "Time Series Plot of 100 Random Numbers", xlab = "Sample Index", ylab = "Value")
hist(Sample100, main = "Histogram of 100 Random Numbers", xlab = "Value", ylab = "Frequency")

# iii)
Sample1000 = Middle_Square(1000)
ts.plot(Sample1000, main = "TS Plot of 1000 Random Numbers", xlab = "Sample Index", ylab = "Value")
hist(Sample1000, main = "Histogram of 1000 Random Numbers", xlab = "Value", ylab = "Frequency")

# iv)
Sample10000 = Middle_Square(10000)
ts.plot(Sample10000, main = "Time Series Plot of 10000 Random Numbers", xlab = "Sample Index", ylab = "Value")
hist(Sample10000, main = "Histogram of 10000 Random Numbers", xlab = "Value", ylab = "Frequency")

Sample15000 = Middle_Square(15000)
ts.plot(Sample15000, main = "Time Series Plot of 15000 Random Numbers", xlab = "Sample Index", ylab = "Value")
hist(Sample15000, main = "Histogram of 15000 Random Numbers", xlab = "Value", ylab = "Frequency")