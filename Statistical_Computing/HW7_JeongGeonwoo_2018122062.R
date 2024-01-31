# STA3127 Statistical Computing
# HW7
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)

##### Q1 #####
### i) ###
my_rmvnorm <- function(n, mean, sigma) {
  # Cholesky decomposition
  L <- chol(sigma) 
  # Box-Muller Transformation
  my_norm <- function(n) { 
    u1 <- runif(n/2); u2 <- runif(n/2)
    z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)
    return(c(z1, z2))}
  # Standard Normal distribution
  z <- matrix(my_norm(n * length(mean)), ncol = length(mean), byrow = TRUE)
  # cov matrix
  mean + z %*%L}

mvnorm_density <- function(x_vec, mean_vec, cov_matrix) {
  k <- length(mean_vec)
  sqrt_det_cov <- sqrt(det(cov_matrix))
  inv_cov_matrix <- solve(cov_matrix)
  diff_vec <- x_vec - mean_vec
  exponent <- -0.5 * t(diff_vec) %*% inv_cov_matrix %*% diff_vec
  return( (1 / ((2 * pi)^(k / 2) * sqrt_det_cov)) * exp(exponent)[1,1] )}

f.density <- function(vec){
  return(mvnorm_density(vec, c(0,0), matrix(c(1,0.5,0.5,1),ncol=2)))}

g.density <- function(vec, v){
  return(mvnorm_density(vec, c(10,10), v*matrix(c(1,0.5,0.5,1),ncol=2)))}

n = 1e6

g.v.0.1 <- my_rmvnorm(n, c(10,10), 0.1 * matrix(c(1, 0.5, 0.5, 1), nrow=2))
g.v.1 <- my_rmvnorm(n, c(10,10), 1 * matrix(c(1, 0.5, 0.5, 1), nrow=2))
g.v.10 <- my_rmvnorm(n, c(10,10), 10 * matrix(c(1, 0.5, 0.5, 1), nrow=2))
g.v.100 <- my_rmvnorm(n, c(10,10), 100 * matrix(c(1, 0.5, 0.5, 1), nrow=2))

indicator <- function(mat){
  # input : n*2 matrix
  # output : n indicator
  return(rowSums(matrix(c((mat[,1] - 10)^2, (mat[,2] - 10)^2), ncol=2))<10)}

weight.v.0.1 <- apply(g.v.0.1, 1, f.density) / apply(g.v.0.1, 1, function(row) g.density(row, 0.1))
weight.v.1 <- apply(g.v.1, 1, f.density) / apply(g.v.1, 1, function(row) g.density(row, 1))
weight.v.10 <- apply(g.v.10, 1, f.density) / apply(g.v.10, 1, function(row) g.density(row, 10))
weight.v.100 <- apply(g.v.100, 1, f.density) / apply(g.v.100, 1, function(row) g.density(row, 100))

cat("v = 0.1 |mean:", mean(indicator(g.v.0.1) * weight.v.0.1), "\n",
    "v = 1 |mean:", mean(indicator(g.v.1) * weight.v.1), "\n",
    "v = 10 |mean:", mean(indicator(g.v.10) * weight.v.10), "\n",
    "v = 100 |mean:", mean(indicator(g.v.100) * weight.v.100))

cat("v = 0.1 |standard error:", sd(indicator(g.v.0.1) * weight.v.0.1) / sqrt(n), "\n",
    "v = 1 |standard error:", sd(indicator(g.v.1) * weight.v.1) / sqrt(n), "\n",
    "v = 10 |standard error:", sd(indicator(g.v.10) * weight.v.10) / sqrt(n), "\n",
    "v = 100 |standard error:", sd(indicator(g.v.100) * weight.v.100) / sqrt(n))

##### Q2 #####
### ii) ###
mu=1; sigma=1; a=1; b=2; K=10; v=1; s=3; t=10; n=1e5

# C(K,t,v)
C <- function(mu, sigma, K, t, v){
  m <- (t * mu - log(K/v)) / sigma / sqrt(t)
  return (v * exp( mu * t + sigma^2 * t / 2)
          + pnorm(sigma * sqrt(t) + m) - K * pnorm(m))}

# ( log(a/v)-t\mu ) / ( sigma * sqrt(s) )
log_normalized_a.v <- (log(a/v)-s*mu) / sigma / sqrt(s)
log_normalized_b.v <- (log(b/v)-s*mu) / sigma / sqrt(s)

# Generating the X from Truncated Normal distribution
q <- (pnorm(log_normalized_b.v) - pnorm(log_normalized_a.v)
      ) * runif(n) + pnorm(log_normalized_a.v)
X_trunc <-  qnorm(q) * sigma * sqrt(s) + s * mu

# Calculate the expected payoff of the option
E.R.Expectation_part = C(mu=mu, sigma=sigma, K=K, t=t-s, v=v*exp(X_trunc))
E.R.CDF_part = pnorm( (log(b/v) - s*mu) / sigma / sqrt(s) ) - pnorm( (log(a/v) - s*mu) / sigma / sqrt(s) )
E.R <- E.R.Expectation_part * E.R.CDF_part

# Check the estimated value and standard error
cat("Expected payoff of the barrier call option:", mean(E.R), "\n",
    "Standard error of payoff of the barrier call option:", sd(E.R)/sqrt(n))

### iii) ###
X <- qnorm(runif(n)) * sqrt(s * sigma^2) + s * mu
Y <- qnorm(runif(n)) * sqrt((t-s) * sigma^2) + (t-s) * mu

# if the condition is not met, payoff becomes zero.
# if the condition is met, payoff becomes v*exp(X+Y)-K
E.R.usual <- ifelse(X <= log(a/v) | X >= log(b/v), 0, v*exp(X+Y)-K)

# Check the estimated value and standard error
cat("Expected payoff of the barrier call option:", mean(E.R.usual), "\n",
    "Standard error of payoff of the barrier call option:", sd(E.R.usual)/sqrt(n))
