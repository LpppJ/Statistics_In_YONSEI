# STA3127 Statistical Computing
# HW5
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
dev.off()
set.seed(2018122062)

##### Q1 #####
### i) ###
# function 1. Define the function to generate Gaussian distribution
# (using Cholesky decomposition and Box-Muller Transformation)
my_mvrnorm <- function(n, mean, sigma) {
  L <- chol(sigma) # Cholesky decomposition
  my_norm <- function(n) { # Box-Muller Transformation
    u1 <- runif(n/2); u2 <- runif(n/2)
    z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)
    return(c(z1, z2))
  }
  # Standard Normal distribution
  z <- matrix(my_norm(n * length(mean)), ncol = length(mean), byrow = TRUE)
  mean + z %*%L # cov matrix
}
# function 2. Define the function to generate Gaussian copula
my_gaussian_copula <- function(n, rho) {
  sigma <- matrix(c(1, rho, rho, 1), nrow=2) # Correlation matrix
  normal_samples <- my_mvrnorm(n, mean=c(0, 0), sigma=sigma) # Generate bivariate normal
  U <- pnorm(normal_samples[,1]); V <- pnorm(normal_samples[,2])
  return(list(U = U, V = V))
}

# Generate samples for rho = -0.5
gaussian.samples.neg <- my_gaussian_copula(1e5, -0.5)
# Histograms for rho = -0.5
par(mfrow=c(1,2))
hist(gaussian.samples.neg$U, xlab="U", breaks=50, freq=FALSE,
     main="Histogram of U (rho = -0.5)")
hist(gaussian.samples.neg$V, xlab="V", breaks=50, freq=FALSE,
     main="Histogram of V (rho = -0.5)")

# Scatter plot for rho = -0.5
par(mfrow=c(1,1))
plot(gaussian.samples.neg$U, gaussian.samples.neg$V,
     main="Scatter plot of U vs V (rho = -0.5)",
     xlab="U", ylab="V", pch=19, cex=0.1)

# Generate samples for rho = 0.5
gaussian.samples.pos <- my_gaussian_copula(1e5, 0.5)
# Histograms for rho = 0.5
par(mfrow=c(1,2))
hist(gaussian.samples.pos$U, xlab="U", breaks=50, freq=FALSE,
     main="Histogram of U (rho = 0.5)")
hist(gaussian.samples.pos$V, xlab="V", breaks=50, freq=FALSE,
     main="Histogram of V (rho = 0.5)")

# Scatter plot for rho = 0.5
par(mfrow=c(1,1))
plot(gaussian.samples.pos$U, gaussian.samples.pos$V,
     main="Scatter plot of U vs V (rho = 0.5)",
     xlab="U", ylab="V", pch=19, cex=0.1)

### ii) ###
# Function to transform
transform_to_exponential <- function(U, V, lambda1, lambda2) {
  X <- -log(1 - U) / lambda1
  Y <- -log(1 - V) / lambda2
  return(list(X = X, Y = Y))
}

# Transform for both rho values
exp.samples.neg <- transform_to_exponential(gaussian.samples.neg$U,
                                            gaussian.samples.neg$V,
                                            2, 3)
exp.samples.pos <- transform_to_exponential(gaussian.samples.pos$U,
                                            gaussian.samples.pos$V, 
                                            2, 3)

# Plot histograms and scatter plots for Exponential marginals
# Histograms for rho = -0.5
par(mfrow=c(1,2))
hist(exp.samples.neg$X, breaks=50, freq=FALSE,
     main="Histogram of X (rho = -0.5)", xlab="X")
hist(exp.samples.neg$Y, xlab="Y", breaks=50, freq=FALSE,
     main="Histogram of Y (rho = -0.5)")

# Scatter plot for rho = -0.5
par(mfrow=c(1,1))
plot(exp.samples.neg$X, exp.samples.neg$Y,
     main="Scatter plot of X vs Y (rho = -0.5)",
     xlab="X", ylab="Y", pch=19, cex=0.1)

# Plot histograms and scatter plots
# Histograms for rho = 0.5
par(mfrow=c(1,2))
hist(exp.samples.pos$X, xlab="X", breaks=50, freq=FALSE,
     main="Histogram of X (rho = 0.5)")
hist(exp.samples.pos$Y, xlab="Y", breaks=50, freq=FALSE,
     main="Histogram of Y (rho = 0.5)")

# Scatter plot for rho = 0.5
par(mfrow=c(1,1))
plot(exp.samples.pos$X, exp.samples.pos$Y,
     main="Scatter plot of X vs Y (rho = 0.5)",
     xlab="X", ylab="Y", pch=19, cex=0.1)

### iii) ###
corr_XY.neg <- cor(exp.samples.neg$X, exp.samples.neg$Y)
corr_XY.pos <- cor(exp.samples.pos$X, exp.samples.pos$Y)
cat("rho : -0.5 ||| sample correlation :", corr_XY.neg, "\n",
    "rho : 0.5 ||| sample correlation :", corr_XY.pos)


##### Q2 #####
n <- 10; r <- 1000
my_rexp <- function(n, lambda){
  U_matrix <- runif(n)
  return (-log(1-U_matrix)/lambda)
}
exp_10 <- my_rexp(n, 1)

### i) ###
index_matrix <- matrix(floor( n * runif(r * n)) +1, nrow = r, ncol = n)
bootstrap_samples <- matrix(exp_10[index_matrix], nrow = r, ncol = n)
bootstrap_var <- apply(bootstrap_samples, 1, var)

hist(bootstrap_var, main = "Bootstrap Sample Variances",
     xlab = "sample variances", xlim = range(c(0, 3)), breaks=50, freq=FALSE, col="red")

### ii) ###
monte_carlo_samples <- matrix(my_rexp(n*r, 1), nrow=r, ncol=n)
monte_carlo_var <- apply(monte_carlo_samples, 1, var)

hist(bootstrap_var, main = "Bootstrap and Monte Carlo Sample Variances",
     xlab = "sample variances", xlim = range(c(0, 3)), breaks=50, freq=FALSE, col="red")
hist(monte_carlo_var, add=TRUE, col="blue", breaks=200, freq=FALSE)
legend("topright", legend = c("Bootstrap", "Monte Carlo"), fill = c("red", "blue"))

### iii) ###
mu4 <- 9
sigma2 <- 1
hist(bootstrap_var, main = "Bootstrap and Monte Carlo Sample Variances",
     xlab = "sample variances", xlim = range(c(0, 3)), breaks=50, freq=FALSE, col="red")
hist(monte_carlo_var, add=TRUE, breaks=200, freq=FALSE, col="blue")
# green line : limiting sampling distribution using mu4=9, sigma2=1
curve(dnorm(x, mean = sigma2, sd=sqrt(abs(mu4-sigma2^2)/n)), col = "green", add = TRUE)
legend("topright", legend = c("Bootstrap", "Monte Carlo","Limiting"), 
       fill = c("red", "blue", "green"))

### iv) ###
hist(bootstrap_var, main = "Bootstrap and Monte Carlo Sample Variances",
     xlab = "sample variances", xlim = range(c(0, 3)), breaks=50, freq=FALSE, col="red")
hist(monte_carlo_var, add=TRUE, col="blue", breaks=200, freq=FALSE)

# green line : limiting sampling distribution using [mu4=9, sigma2=1]
curve(dnorm(x, mean = sigma2, sd=sqrt(abs(mu4-sigma2^2)/n)), col = "green", add = TRUE)

# orange line : approximate sampling distribution using [bootstrap sampling distribution]
mu4 <- mean((bootstrap_samples-mean(bootstrap_samples))^4)
sigma2 <- var(as.vector(bootstrap_samples))
curve(dnorm(x, mean = sigma2, sd=sqrt(abs(mu4-sigma2^2)/n)), col = "orange", add = TRUE)

# legend
legend("topright", legend = c("Bootstrap", "Monte Carlo","Limiting", "Approx_bootstrap"),
       fill = c("red", "blue", "green", "orange"))











### vi) ###
n <- 1e5; r <- 1000
my_rexp <- function(n, lambda){
  U_matrix <- runif(n)
  return (-log(1-U_matrix)/lambda)
}
exp_1e5 <- my_rexp(n, 1)

### vi-i) ###
index_matrix <- matrix(floor( n * runif(r * n)) +1, nrow = r, ncol = n)
bootstrap_samples <- matrix(exp_1e5[index_matrix], nrow = r, ncol = n)
bootstrap_var <- apply(bootstrap_samples, 1, var)

hist(bootstrap_var, main = "Bootstrap Sample Variances",
     xlab = "sample variances", xlim = range(c(0.95, 1.05)), breaks=50, freq=FALSE, col="red")

### vi-ii) ###
monte_carlo_samples <- matrix(my_rexp(n*r, 1), nrow=r, ncol=n)
monte_carlo_var <- apply(monte_carlo_samples, 1, var)

hist(bootstrap_var, main = "Bootstrap and Monte Carlo Sample Variances",
     xlab = "sample variances", xlim = range(c(0.95, 1.05)), breaks=50, freq=FALSE, col="red")
hist(monte_carlo_var, add=TRUE, col="blue", breaks=80, freq=FALSE)
legend("topright", legend = c("Bootstrap", "Monte Carlo"), fill = c("red", "blue"))

### vi-iii) ###
mu4 <- 9
sigma2 <- 1
hist(bootstrap_var, main = "Bootstrap and Monte Carlo Sample Variances",
     xlab = "sample variances", xlim = range(c(0.95, 1.05)), breaks=50, freq=FALSE, col="red")
hist(monte_carlo_var, add=TRUE, col="blue", breaks=80, freq=FALSE)
# green line : limiting sampling distribution using mu4=9, sigma2=1
curve(dnorm(x, mean = sigma2, sd=sqrt(abs(mu4-sigma2^2)/n)), col = "green", add = TRUE)
legend("topright", legend = c("Bootstrap", "Monte Carlo","Limiting"), 
       fill = c("red", "blue", "green"))

### vi-iv) ###
hist(bootstrap_var, main = "Bootstrap and Monte Carlo Sample Variances",
     xlab = "sample variances", xlim = range(c(0.95, 1.05)), breaks=50, freq=FALSE, col="red")
hist(monte_carlo_var, add=TRUE, col="blue", breaks=80, freq=FALSE)

# green line : limiting sampling distribution using [mu4=9, sigma2=1]
curve(dnorm(x, mean = sigma2, sd=sqrt(abs(mu4-sigma2^2)/n)), col = "green", add = TRUE)

# orange line : approximate sampling distribution using [bootstrap sampling distribution]
mu4 <- mean((bootstrap_samples-mean(bootstrap_samples))^4)
sigma2 <- var(as.vector(bootstrap_samples))
curve(dnorm(x, mean = sigma2, sd=sqrt(abs(mu4-sigma2^2)/n)), col = "orange", add = TRUE)

# legend
legend("topright", legend = c("Bootstrap", "Monte Carlo","Limiting", "Approx_bootstrap"),
       fill = c("red", "blue", "green", "orange"))
