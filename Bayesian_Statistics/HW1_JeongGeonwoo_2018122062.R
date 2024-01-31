##### Problem 1 #####
# Setting and sampling
library(coda)
set.seed(2018122062)
mean <- 189/20
sd <- sqrt(3/4)
samples <- rnorm(1000000, mean = mean, sd = sd)

# Posterior mean
cat("Posterior mean :",
    mean(samples)
    )

# Highest Posterior Density (HPD) Interval
mcmc_samples <- as.mcmc(samples)
hpd_interval <- HPDinterval(mcmc_samples, prob = 0.95)
cat("95% highest posterior density interval
    [",hpd_interval,"]"
    )

# Symmetrical Density Interval
sorted_samples <- sort(samples)
lower_percentile <- 0.025
upper_percentile <- 0.975
lower_index <- floor(length(sorted_samples) * lower_percentile) + 1
upper_index <- ceiling(length(sorted_samples) * upper_percentile)
confidence_interval <- c(sorted_samples[lower_index],
                         sorted_samples[upper_index])
cat("95% symmetrical density interval
    [", confidence_interval[1], ", ", confidence_interval[2],"]"
    )

##### Problem 2 #####
# Setting and sampling
set.seed(2018122062)
shape <- 11/2
scale <- 21
gamma_samples <- rgamma(1000000, shape=shape, scale=1/scale)
ig_samples <- 1 / gamma_samples

# Posterior mean
cat("Posterior mean :",
    mean(ig_samples)
    )

# Highest Posterior Density (HPD) Interval
mcmc_ig_samples <- as.mcmc(ig_samples)
hpd_interval <- HPDinterval(mcmc_ig_samples, prob = 0.95)
cat("95% highest posterior density interval
    [",hpd_interval,"]"
    )

# Symmetrical Density Interval
sorted_ig_samples <- sort(ig_samples)
lower_percentile <- 0.025
upper_percentile <- 0.975
lower_index <- floor(length(sorted_ig_samples) * lower_percentile) + 1
upper_index <- ceiling(length(sorted_ig_samples) * upper_percentile)
confidence_interval <- c(sorted_ig_samples[lower_index],
                         sorted_ig_samples[upper_index])
cat("95% symmetrical density interval
    [", confidence_interval[1], ", ", confidence_interval[2],"]"
    )
