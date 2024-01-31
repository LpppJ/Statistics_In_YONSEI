# STA3127 Statistical Computing
# HW6
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)

##### Q1 #####
### iii) ###
r=10; p=0.3

X <- function(u.vec, p=0.3){
  r <- length(u.vec)
  x <- sum( ceiling( log(u.vec) / log(p) ) ) - r
  return (x)}

b <- function(X){return (log(exp(X)+1))}

# I will start with 1st sample
n.b <- 1
Bn.bar <- b(X(runif(r)))
Bn.s2 <- 0
while (TRUE){
  # 2rd sample is always needed.
  # Why? 1 : CLT is established at n > 30
  # Why? 2 : If 2rd sample is not needed, sample variance = 0
  #          then the conditional statement below cannot be calculated.
  next_b <- b(X(runif(r)))
  # successive computation for computation speed
  new.Bn.bar <- Bn.bar + (next_b - Bn.bar) / (n.b+1)
  Bn.s2 <- (1-1/n.b) * Bn.s2 + (n.b+1) * (new.Bn.bar - Bn.bar)^2
  Bn.bar <- new.Bn.bar
  n.b <- n.b + 1
  # 3,4,...th samples is conditionally needed.
  # If sample size is enough, stop sampling
  if (0.01 / sqrt(Bn.s2/(n.b)) > 1.959964){break}}

cat("n of B", n.b, "\n", "Bn_bar", Bn.bar)

### iv) ###
V <- function(u.vec, p=0.3){
  return ( ( b(X(u.vec)) + b(X(1-u.vec)) ) / 2 )
}

# I will start with 1st sample
n.v<- 1
Vn.bar <- V(runif(r))
Vn.s2 <- 0
while (TRUE){
  # 2rd sample is always needed.
  # Why? 1 : CLT is established at n > 30
  # Why? 2 : If 2rd sample is not needed, sample variance = 0
  #          then the conditional statement below cannot be calculated.
  next_V <- V(runif(r))
  # successive computation for computation speed
  new.Vn.bar <- Vn.bar + (next_V - Vn.bar) / (n.v+1)
  Vn.s2 <- (1-1/n.v) * Vn.s2 + (n.v+1) * (new.Vn.bar - Vn.bar)^2
  Vn.bar <- new.Vn.bar
  n.v <- n.v + 1
  # 3,4,...th samples is conditionally needed.
  # If sample size is enough, stop sampling
  if (0.01 / sqrt(Vn.s2/(n.v)) > 1.959964){break}}

cat("n of B", n.v, "\n", "Bn_bar", Vn.bar)

### v) ### 
W <- function(u.vec, p=0.3){
  return ( b(X(u.vec)) - (sum(log(u.vec))+r)/log(p) )
}

# I will start with 1st sample
n.w<- 1
Wn.bar <- W(runif(r))
Wn.s2 <- 0
while (TRUE){
  # 2rd sample is always needed.
  # Why? 1 : CLT is established at n > 30
  # Why? 2 : If 2rd sample is not needed, sample variance = 0
  #          then the conditional statement below cannot be calculated.
  next_W <- W(runif(r))
  # successive computation for computation speed
  new.Wn.bar <- Wn.bar + (next_W - Wn.bar) / (n.w+1)
  Wn.s2 <- (1-1/n.w) * Wn.s2 + (n.w+1) * (new.Wn.bar - Wn.bar)^2
  Wn.bar <- new.Wn.bar
  n.w <- n.w + 1
  # 3,4,...th samples is conditionally needed.
  # If sample size is enough, stop sampling
  if (0.01 / sqrt(Wn.s2/(n.w)) > 1.959964){break}}

cat("n of B", n.w, "\n", "Bn_bar", Wn.bar)
num_simulation<-100; n<-1e4; mean<-1; sd=2
##### Q2 #####
### iv) ###
estimate_usual <- function(num_simulation, n, mean, sd){
  # I will define "usual" function which is
  # estimator E(log(exp(X)+1)) one time
  usual <- function(n, mean, sd){
    normal_samples <- mean + sd * qnorm(runif(n))
    estimator <- mean(log(exp(normal_samples) + 1))
    return (estimator)
  }
  # I repeated usual function "num_simulation" times and get variance
  estimators <- replicate(num_simulation, usual(n, mean, sd))
  result <- list(mean = mean(estimators), var = var(estimators))
  return (result)
}

estimate_stratified <- function(num_simulation, n, mean, sd) {
  # I divided (0,1) into 1e4 intervals
  #   and generated one number from each interval using runif(interval)
  break_points <- seq(0, 1, length.out = n + 1)
  estimators <- replicate(num_simulation, {
    unif_samples <- runif(n, min = head(break_points, -1), max = tail(break_points, -1))
    normal_samples <- mean + sd * qnorm(unif_samples)
    # normal_samples is the samples from N(-1,2^2).
    mean(log(exp(normal_samples) + 1))
  })
  
  result <- list(mean = mean(estimators), var = var(estimators))
  return (result)
}

ptm <- proc.time()
result.usual <- estimate_usual(1e5, 1e4, -1, 2)
elapsed.time.usual <- proc.time() - ptm

ptm <- proc.time()
result.stratified <- estimate_stratified(1e5, 1e4, -1, 2)
elapsed.time.stratified <- proc.time() - ptm

cat("The variance of usual average estimator|",
    "mean:", result.usual$mean, "variance:", result.usual$var, "\n",
    "The variance of stratified estimator|",
    "mean:", result.stratified$mean, "variance:", result.stratified$var)

elapsed.time.usual; elapsed.time.stratified