# STA3127 Statistical Computing
# HW8
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)
library(coda); library(scatterplot3d)
##### Q1 #####
### i) ###
g <- function(x, y) {
  exp(-(x + 1)^2 - y^2) + exp(-150*(x^2 - y)^2 - 150*(x - y^2)^2)}

x <- seq(-3, 2, length.out = 100);y <- seq(-2, 2, length.out = 100)
z <- outer(x, y, g)
par(mfrow=c(1,1))
contour(x, y, z, xlab = "X", ylab = "Y",
        main="The contour plot of density g")

# Additional visualization for understanding the density of g
par(mfrow=c(1,2))
x <- seq(-3, 2, length.out = 100); y <- seq(-2, 2, length.out = 100)
z <- outer(x, y, g)
res1 <- persp(x, y, z, theta = 30, phi = 20, expand = 0.5,
              col = "green", xlab = "X", ylab = "Y", zlab = "Z",
              main = "The density plot of g")
res2 <- persp(x, y, z, theta = 70, phi = 20, expand = 0.5,
              col = "green", xlab = "X", ylab = "Y", zlab = "Z",
              main = "The density plot of g")

par(mfrow=c(1,2))
x <- seq(-0.2, 0.2, length.out = 100); y <- seq(-0.2, 0.2, length.out = 100)
z <- outer(x, y, g)
contour(x, y, z, xlab = "X", ylab = "Y",
        main="Density of g : x,y=(0,0)")
x <- seq(0.9, 1.1, length.out = 100); y <- seq(0.9, 1.1, length.out = 100)
z <- outer(x, y, g)
contour(x, y, z, xlab = "X", ylab = "Y",
        main="Density of g : x,y=(1,1)")

### ii) ###
# Define the function to generate Gaussian distribution
# (using Cholesky decomposition and Box-Muller Transformation)
my_rmvnorm <- function(n, mean, sigma) {
  L <- chol(sigma) # Cholesky decomposition
  my_norm <- function(n) { # Box-Muller Transformation
    u1 <- runif(n/2); u2 <- runif(n/2)
    z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)
    return(c(z1, z2))}
  # Standard Normal distribution
  z <- matrix(my_norm(n * length(mean)), ncol = length(mean), byrow = TRUE)
  # cov matrix
  mean + z %*%L}

# Implement the Metropolis-Hastings algorithm
MCMC.MH <- function(n, init, sigma, rho) {
  x <- rep(NA, n); y <- rep(NA, n)
  x[1] <- init[1]; y[1] <- init[2]
  
  for (i in 2:n) {
    current <- c(x[i-1], y[i-1])
    proposal <- my_rmvnorm(1, mean = current,
                           sigma = matrix(c(1, rho, rho, 1) * sigma^2, ncol = 2))
    a <- g(proposal[1], proposal[2]) / g(current[1], current[2])
    
    if (runif(1) < a) {x[i] <- proposal[1]; y[i] <- proposal[2]}
    else {x[i] <- current[1]; y[i] <- current[2]}}
  return(list(x = x, y = y))}

n <- 1e6
samples.1 <- MCMC.MH(n = n, init = c(-10, -10), sigma = 0.1, rho = 0)
par(mfrow=c(1,1))
plot(samples.1$x, samples.1$y, main = "Scatter plot of x and y", pch = 20, cex = 0.1)

par(mfrow = c(1,2))
hist(samples.1$x, main = "Histogram of x", freq=FALSE, breaks=100)
hist(samples.1$y, main = "Histogram of y", freq=FALSE, breaks=100)
ts.plot(samples.1$x, main = "Trace plot of x")
ts.plot(samples.1$y, main = "Trace plot of y")
acf(samples.1$x, main = "Autocorrelation of x")
acf(samples.1$y, main = "Autocorrelation of y")

### iii) ###
samples.2 <- MCMC.MH(n = n, init = c(-10, -10), sigma = 1, rho = 0.8)
par(mfrow=c(1,1))
plot(samples.2$x, samples.2$y, main = "Scatter plot of x and y", pch = 20, cex = 0.1)

par(mfrow = c(1,2))
hist(samples.2$x, main = "Histogram of x", freq=FALSE, breaks=100)
hist(samples.2$y, main = "Histogram of y", freq=FALSE, breaks=100)
ts.plot(samples.2$x, main = "Trace plot of x")
ts.plot(samples.2$y, main = "Trace plot of y")
acf(samples.2$x, main = "Autocorrelation of x")
acf(samples.2$y, main = "Autocorrelation of y")

### iv) ###
samples.3 <- MCMC.MH(n = n, init = c(-10, -10), sigma = 1, rho = -0.8)
par(mfrow=c(1,1))
plot(samples.3$x, samples.3$y, main = "Scatter plot of x and y", pch = 20, cex = 0.1)

par(mfrow = c(1,2))
hist(samples.3$x, main = "Histogram of x", freq=FALSE, breaks=100)
hist(samples.3$y, main = "Histogram of y", freq=FALSE, breaks=100)
ts.plot(samples.3$x, main = "Trace plot of x")
ts.plot(samples.3$y, main = "Trace plot of y")
acf(samples.3$x, main = "Autocorrelation of x")
acf(samples.3$y, main = "Autocorrelation of y")

### vii) ###

samples.4 <- MCMC.MH(n = n, init = c(-10, -10), sigma = 4, rho = 0.8)
par(mfrow=c(1,1))
plot(samples.4$x, samples.4$y, main = "Scatter plot of x and y", pch = 20, cex = 0.1)

for (i in 1:10){points(samples.4$x[i], samples.4$y[i], col = i, pch = 20, cex = 1)}

par(mfrow = c(1,2))
hist(samples.4$x, main = "Histogram of x", freq=FALSE, breaks=100)
hist(samples.4$y, main = "Histogram of y", freq=FALSE, breaks=100)
ts.plot(samples.4$x, main = "Trace plot of x")
ts.plot(samples.4$y, main = "Trace plot of y")
acf(samples.4$x, main = "Autocorrelation of x")
acf(samples.4$y, main = "Autocorrelation of y")

##### Q2 #####
### i) ###
rm(list=ls())

# Random generating From truncated exponential distribution
r.trunc.exp <- function(rate, lower, upper) {
  # CDF of exponential distribution at lower and upper
  F_lower <- 1 - exp(-rate * lower)
  F_upper <- 1 - exp(-rate * upper)
  
  # runif from (F_lower, F_upper)
  u <- runif(1) * (F_upper - F_lower) + F_lower
  
  # Truncated sample by inverse CDF of exponential distribution
  x <- -log(1 - u) / rate
  
  return(x)
}

# arbitrary initial value
X <- rep(0.5, 3)
# sample size (iteration)
n <- 1e6
# burn in size
burn_in <- 1e2

# sample space
samples <- matrix(NA, nrow = n, ncol = 3)

# Gibbs sampling
for (i in 1:n) {
  a1 <- sqrt(max(0, 1 - 2*X[2]^2 - 3*X[3]^2))
  c1 <- sqrt(2 - 2*X[2]^2 - 3*X[3]^2) - a1
  X[1] <- r.trunc.exp(1.5, 0, c1) + a1
  
  a2 <- sqrt(max(0, (1 - X[1]^2 - 3*X[3]^2)/2))
  c2 <- sqrt((2 - X[1]^2 - 3*X[3]^2)/2) - a2
  X[2] <- r.trunc.exp(1.5, 0, c2) + a2
  
  a3 <- sqrt(max(0, (1 - X[1]^2 - 2*X[2]^2)/3))
  c3 <- sqrt((2 - X[1]^2 - 2*X[2]^2)/3) - a3
  X[3] <- r.trunc.exp(1.5, 0, c3) + a3
  samples[i,] <- X
}

# burn-in
samples <- samples[(burn_in + 1):n, ]

# scatter plot
par(mfrow=c(1,2))
s3d <- scatterplot3d(samples[,1], samples[,2], samples[,3],
              main="3D Scatter Plot : X1, X2, X3",
              xlab="X1", ylab="X2", zlab="X3",
              color="blue", pch=20, angle=45)
s3d$points3d(samples[1:10,1], samples[1:10,2], samples[1:10,3],
             col="red", pch=19)
s3d2 <- scatterplot3d(samples[,1], samples[,2], samples[,3],
                     main="3D Scatter Plot : X1, X2, X3",
                     xlab="X1", ylab="X2", zlab="X3",
                     color="blue", pch=20, angle=225)
s3d2$points3d(samples[1:10,1], samples[1:10,2], samples[1:10,3],
             col="red", pch=19)

# mean vector and covariance matrix
mean_vector <- apply(samples, 2, mean)
covariance_matrix <- cov(samples)

print(mean_vector)
print(covariance_matrix)