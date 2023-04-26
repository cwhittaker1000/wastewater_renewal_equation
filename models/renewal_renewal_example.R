# Loading required libraries
duration <- 350
initial_infections <- 10
SI_cutoff <- 200
N0 <- 10
N2 <- duration
infections <- vector(mode = "numeric", length = duration)
deaths <- vector(mode = "numeric", length = duration)
cumm_sum <- vector(mode = "numeric", length = duration)
infection_SI <- dgamma(x = 1:N2, shape = 2, scale = 1 / 0.1)
infection_SI_rev <- rev(infection_SI)
death_SI <- dgamma(x = 1:N2, shape = 2, scale = 1 / 0.1)
death_SI_rev <- rev(death_SI)
Death_Fitting_Start <- 10
Rt <- c(rep(2.5, 150), rep(1.2, 100), rep(0.9, 50), rep(0.5, 50))
Rt_adj <- vector(mode = "numeric", length = duration)
population <- 10000
infections[1:N0] <- initial_infections

for (i in N0:duration) {
  if ((i - SI_cutoff - 1) < 1) {
    convolution <- sum(infections[1:(i - 1)] * tail(infection_SI_rev, i - 1)) 
  } else {
    convolution <- sum(infections[(i - SI_cutoff):(i - 1)] * tail(infection_SI_rev, SI_cutoff)) 
  }
  cumm_sum[i] <- cumm_sum[i - 1] + infections[i - 1]
  Rt_adj[i] <- Rt[i] # ((pop - cumm_sum[i]) / pop) * Rt[i]
  infections[i] <- Rt_adj[i] * convolution
}

ifr <- 0.01
for (i in 2:duration) {
  if ((i - SI_cutoff - 1) < 1) {
    deaths[i] <- ifr * sum(infections[1:(i - 1)] * tail(death_SI_rev, i - 1)) 
  } else {
    deaths[i] <- ifr * sum(infections[(i - SI_cutoff):(i - 1)] * tail(death_SI_rev, SI_cutoff)) 
  }
}

# Plotting outputs and generating data for model fitting
par(mfrow = c(1, 3))
plot(Rt, type = "l", ylim = c(0, max(Rt)))
plot(infections, type = "l")
set.seed(11)
noisy_data <- rnbinom(n = duration, mu = deaths, size = 50) # low size = more overdispersion
matplot(deaths, xlab = "Time", ylab = "Number of deaths", type = "l", col = "blue", lty = 1, ylim = c(0, max(noisy_data)))
matpoints(noisy_data, col = "blue", pch = 20)

# Generating data list
week_index <- c(rep(1, sum(Rt == unique(Rt)[1])),
                rep(2, sum(Rt == unique(Rt)[2])),
                rep(3, sum(Rt == unique(Rt)[3])),
                rep(4, sum(Rt == unique(Rt)[3])))
data_list <- list(N0 = N0, 
                  N2 = N2,
                  ifr_mean = 0.02,
                  ifr_sd = 0.02,
                  deaths = noisy_data,
                  Death_Fitting_Start = 10,
                  week_index = week_index, 
                  W = max(week_index),
                  SI_cutoff = 100,
                  SI = infection_SI,
                  death_delay = death_SI)

# try to work out what's going on with the weird infections values (negatives)
mod <- cmdstanr::cmdstan_model(stan_file = "models/base_renewal_equation_NoPopImm.stan")
fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 100,
  list(
    list(mu = 2,
         kappa = 0.5,
         initial_infections = 2,
         tau = 5,
         phi = 5))
)
draws_df <- fit$draws(format = "df")
fit$summary()
par(mfrow = c(3, 2))
hist(draws_df$mu)
hist(draws_df$kappa)
hist(draws_df$ifr)
hist(draws_df$phi)
hist(draws_df$initial_infections)
hist(draws_df$tau)

# Infections Over Time
par(mfrow = c(1, 3))
infections_output <- draws_df[, grep("infections0*", colnames(draws_df))[1:duration]]
infections <- unname(colMeans(infections_output))
infections_lower <- unname(apply(infections_output, 2, quantile, 0.05))
infections_upper <- unname(apply(infections_output, 2, quantile, 0.95))

matplot(infections, xlab = "Time", ylab = "Number of infections", pch = 20, col = "black", 
        lty = 1, ylim = c(0, max(infections_upper[1:duration])))
matlines(infections[1:duration], col = "red", pch = 20)
matlines(infections_lower[1:duration], col = "red", pch = 20)
matlines(infections_upper[1:duration], col = "red", pch = 20)

# deaths Over Time
deaths_output <- draws_df[, grep("deaths0*", colnames(draws_df))[1:duration]]
deaths <- unname(colMeans(deaths_output))
deaths_lower <- unname(apply(deaths_output, 2, quantile, 0.05))
deaths_upper <- unname(apply(deaths_output, 2, quantile, 0.95))

matplot(noisy_data[1:duration], xlab = "Time", ylab = "Number of deaths", col = "black", pch = 20,
        ylim = c(0, max(c(noisy_data[1:duration], deaths[1:duration]))))
matlines(deaths[1:duration], col = "red", pch = 20)
matlines(deaths_lower[1:duration], col = "red", pch = 20)
matlines(deaths_upper[1:duration], col = "red", pch = 20)

# Rt over time
Rt <- draws_df[grep("Rt\\[", colnames(draws_df))]
plot(apply(Rt, 2, quantile, 0.95)[1:duration], type = "l", ylab = "Reff", ylim = c(0, max(apply(Rt, 2, quantile, 0.95))))
lines(apply(Rt, 2, quantile, 0.05)[1:duration], type = "l", ylab = "Reff")
lines(colMeans(Rt)[1:duration], type = "l", ylab = "Reff")
